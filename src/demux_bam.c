#include "../include/common.h"
#include "../include/utils.h"
#include <htslib/sam.h>

#if GLIB_CHECK_VERSION(2,68,0)
#define MEMDUP g_memdup2
#else
#define MEMDUP g_memdup
#endif

/*
 * demux_bam.c – Read STARsolo BAM, deduplicate (CB,UB,gene) and output MatrixMarket.
 * Phase-1: single consumer ( --threads 1) .
 * Phase-2: optional shard workers for --threads >1 .
 */

#define DEFAULT_CB_TAG "CB"
#define DEFAULT_UB_TAG "UB"
#define DEFAULT_GENE_TAG "GX"
#define FALLBACK_GENE_TAG "GE"

/* ---------------- Intern tables ---------------- */

typedef struct {
    GHashTable *str_to_id;   /* char*  -> guint32 */
    GPtrArray  *id_to_str;   /* index -> char*   */
} intern_table;

static inline guint32 intern_get(intern_table *t, const char *s) {
    gpointer val = g_hash_table_lookup(t->str_to_id, s);
    if (val) return GPOINTER_TO_UINT(val);
    guint32 new_id = t->id_to_str->len;
    g_hash_table_insert(t->str_to_id, g_strdup(s), GUINT_TO_POINTER(new_id));
    g_ptr_array_add(t->id_to_str, g_strdup(s));
    return new_id;
}

static inline const char* id_to_str(intern_table *t, guint32 id) {
    return (const char*)g_ptr_array_index(t->id_to_str, id);
}

static intern_table* intern_table_new(void) {
    intern_table *t = g_new0(intern_table,1);
    t->str_to_id = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
    t->id_to_str = g_ptr_array_new_with_free_func(g_free);
    return t;
}

/* ---------------- Sample map ---------------- */

typedef struct {
    GHashTable *cb_to_sample;   /* char* CB -> guint8 sample_index */
    GHashTable *sample_to_index;/* char* sample_key -> guint8 */
    guint8 next_index;          /* next available (starts at 1; 0 = unassigned) */
} sample_map;

static sample_map* sample_map_new(void) {
    sample_map *m = g_new0(sample_map,1);
    m->cb_to_sample = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
    m->sample_to_index = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
    m->next_index = 1;
    return m;
}

static guint8 sample_map_get_index(sample_map *m, const char *cb) {
    gpointer v = g_hash_table_lookup(m->cb_to_sample, cb);
    if (!v) return 0; /* unassigned */
    return GPOINTER_TO_UINT(v);
}

static void load_cell_sample_map(sample_map *m, const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) { perror("cell_sample_map"); exit(EXIT_FAILURE); }
    char line[LINE_LENGTH];
    while (fgets(line, sizeof(line), fp)) {
        char *p = strchr(line,'\n'); if (p) *p = 0;
        if (!line[0]) continue;
        char *a = strtok(line, "\t");
        char *b = strtok(NULL, "\t");
        char *c = strtok(NULL, "\t");
        const char *cb=NULL,*sample_key=NULL;
        if (c) {   /* 3-column */
            cb = b; sample_key = c; /* lib is a, ignoring for now */
        } else if (a && b) { /* 2-column */
            cb = a; sample_key = b;
        }
        if (!cb || !sample_key) continue;
        gpointer idx_ptr = g_hash_table_lookup(m->sample_to_index, sample_key);
        guint8 idx;
        if (!idx_ptr) {
            idx = m->next_index++;
            if (idx == 0) { fprintf(stderr,"Too many samples (>255)\n"); exit(EXIT_FAILURE);} 
            g_hash_table_insert(m->sample_to_index, g_strdup(sample_key), GUINT_TO_POINTER(idx));
        } else idx = GPOINTER_TO_UINT(idx_ptr);
        g_hash_table_insert(m->cb_to_sample, g_strdup(cb), GUINT_TO_POINTER(idx));
    }
    fclose(fp);
    fprintf(stderr,"Loaded %u sample assignments\n", g_hash_table_size(m->cb_to_sample));
}

/* ---------------- 2-bit UMI pack ---------------- */
static inline uint8_t base2bits(char c) {
    switch(c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4; /* invalid */
    }
}

static uint64_t pack_umi_2bit(const char *s, int *ok) {
    uint64_t out = 0; int len=0; *ok=1;
    for (; s[len] && len < 30; ++len) {
        uint8_t bits = base2bits(s[len]);
        if (bits > 3) { *ok=0; return 0;} /* invalid */
        out |= ((uint64_t)bits) << (len*2);
    }
    out |= ((uint64_t)len) << 60; /* store length high bits (up to 15) */
    return out;
}

/* ---------------- TripKey hashing ---------------- */

typedef struct {
    uint64_t cbk;
    uint32_t gene_idx;
    uint64_t umi64;
} TripKey;

static guint tripkey_hash(gconstpointer v) {
    const TripKey *k = v;
    uint64_t h = k->cbk ^ ((uint64_t)k->gene_idx<<32) ^ k->umi64;
    h ^= h >> 33; h *= 0xff51afd7ed558ccdULL; h ^= h >> 33; h *= 0xc4ceb9fe1a85ec53ULL; h ^= h >> 33;
    return (guint)h;
}

static gboolean tripkey_equal(gconstpointer a, gconstpointer b) {
    const TripKey *x = a; const TripKey *y = b;
    return x->cbk==y->cbk && x->gene_idx==y->gene_idx && x->umi64==y->umi64;
}

/* Key for counts (cb_id,gene_idx) packed into uint64 */
static inline uint64_t pair64(uint32_t a, uint32_t b) {
    return ((uint64_t)a<<32)|b;
}

/* ---------------- Shard data ---------------- */

typedef struct {
    GHashTable *dedup_set;   /* TripKey* -> GINT_TO_POINTER(1) */
    GHashTable *counts_map;  /* uint64_t* -> uint32 count (GUINT_TO_POINTER) */
} shard_data;

static void tripkey_free(gpointer data){ g_slice_free(TripKey, data); }

static shard_data* shard_data_new(void) {
    shard_data *s = g_new0(shard_data,1);
    s->dedup_set = g_hash_table_new_full(tripkey_hash, tripkey_equal, tripkey_free, NULL);
    s->counts_map = g_hash_table_new_full(g_int64_hash, g_int64_equal, g_free, NULL);
    return s;
}

/* ---------------- Seen reads table ---------------- */
static guint64 hash_qname(const char *qname, int read12) {
    /* simple 64-bit hash (FNV-1a) */
    uint64_t h=14695981039346656037ULL;
    for (const unsigned char *p=(const unsigned char*)qname; *p; ++p) { h ^= *p; h *= 1099511628211ULL; }
    h ^= read12;
    return h;
}

/* ---------------- CLI config ---------------- */

typedef struct {
    const char *bam_path;
    const char *outdir;
    const char *sample_map_path;
    const char *cb_tag;
    const char *ub_tag;
    const char *gene_tag;
    int hts_threads;
    int consumer_threads; /* shards */
    int min_mapq;
    int primary_only;
    int skip_dup;
    long long max_records;
    int debug;
} cfg_t;

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s --bam FILE --outdir DIR [options]\n"
        "Options:\n"
        "  --cell_sample_map FILE   Optional CB→sample mapping TSV\n"
        "  --cb_tag STR             Default CB\n"
        "  --ub_tag STR             Default UB\n"
        "  --gene_tag STR           Default GX (fallback GE)\n"
        "  --hts_threads N          BGZF threads (default 2)\n"
        "  --threads N              Consumer shards (default 1)\n"
        "  --min_mapq N             Minimum MAPQ (default 0)\n"
        "  --no_primary_filter      Keep secondary/supplementary\n"
        "  --keep_dup               Keep duplicate reads\n"
        "  --max_records N          Process first N records only\n"
        "  -v                       Debug\n",
        prog);
}

static void parse_cli(cfg_t *cfg, int argc, char **argv) {
    static struct option long_opts[] = {
        {"bam", required_argument, 0, 0},
        {"outdir", required_argument,0,0},
        {"cell_sample_map", required_argument,0,0},
        {"cb_tag", required_argument,0,0},
        {"ub_tag", required_argument,0,0},
        {"gene_tag", required_argument,0,0},
        {"hts_threads", required_argument,0,0},
        {"threads", required_argument,0,0},
        {"min_mapq", required_argument,0,0},
        {"no_primary_filter", no_argument,0,0},
        {"keep_dup", no_argument,0,0},
        {"max_records", required_argument,0,0},
        {"debug", no_argument,0,'v'},
        {0,0,0,0}
    };
    int opt_idx=0;
    while (1) {
        int c = getopt_long(argc, argv, "v", long_opts, &opt_idx);
        if (c==-1) break;
        if (c=='v') { cfg->debug=1; continue; }
        if (c!=0) continue;
        const char *opt_name = long_opts[opt_idx].name;
        if (!strcmp(opt_name, "bam")) cfg->bam_path = optarg;
        else if (!strcmp(opt_name, "outdir")) cfg->outdir = optarg;
        else if (!strcmp(opt_name, "cell_sample_map")) cfg->sample_map_path = optarg;
        else if (!strcmp(opt_name, "cb_tag")) cfg->cb_tag = optarg;
        else if (!strcmp(opt_name, "ub_tag")) cfg->ub_tag = optarg;
        else if (!strcmp(opt_name, "gene_tag")) cfg->gene_tag = optarg;
        else if (!strcmp(opt_name, "hts_threads")) cfg->hts_threads = atoi(optarg);
        else if (!strcmp(opt_name, "threads")) cfg->consumer_threads = atoi(optarg);
        else if (!strcmp(opt_name, "min_mapq")) cfg->min_mapq = atoi(optarg);
        else if (!strcmp(opt_name, "no_primary_filter")) cfg->primary_only = 0;
        else if (!strcmp(opt_name, "keep_dup")) cfg->skip_dup = 0;
        else if (!strcmp(opt_name, "max_records")) cfg->max_records = atoll(optarg);
    }
    if (!cfg->bam_path || !cfg->outdir) { usage(argv[0]); exit(EXIT_FAILURE);}    
    if (!cfg->cb_tag) cfg->cb_tag = DEFAULT_CB_TAG;
    if (!cfg->ub_tag) cfg->ub_tag = DEFAULT_UB_TAG;
    if (!cfg->gene_tag) cfg->gene_tag = DEFAULT_GENE_TAG;
    if (cfg->hts_threads<=0) cfg->hts_threads=2;
    if (cfg->consumer_threads<=0) cfg->consumer_threads=1;
}

/* ---------------- Worker function (single-thread) ---------------- */

static void process_bam_single(cfg_t *cfg, intern_table *cb_tab, intern_table *gene_tab,
                               sample_map *smap) {
    samFile *in = sam_open(cfg->bam_path, "r");
    if (!in) { perror("sam_open"); exit(EXIT_FAILURE);}    
    bam_hdr_t *hdr = sam_hdr_read(in);
    hts_set_threads(in, cfg->hts_threads);

    shard_data *shd = shard_data_new();
    GHashTable *seen_reads = g_hash_table_new(g_direct_hash, g_direct_equal);

    bam1_t *b = bam_init1();
    long long recs=0, used=0;
    while (sam_read1(in, hdr, b) >= 0) {
        if (cfg->max_records && recs>=cfg->max_records) break;
        recs++;
        uint16_t flag = b->core.flag;
        if (cfg->primary_only && (flag & (BAM_FSECONDARY|BAM_FSUPPLEMENTARY))) continue;
        if (cfg->skip_dup && (flag & BAM_FDUP)) continue;
        if (b->core.qual < cfg->min_mapq) continue;

        /* Tags */
        const uint8_t *cb_tag = bam_aux_get(b, cfg->cb_tag);
        const uint8_t *ub_tag = bam_aux_get(b, cfg->ub_tag);
        const uint8_t *gx_tag = bam_aux_get(b, cfg->gene_tag);
        if (!gx_tag) gx_tag = bam_aux_get(b, FALLBACK_GENE_TAG);
        if (!cb_tag || !ub_tag || !gx_tag) continue;
        const char *cb = bam_aux2Z(cb_tag);
        const char *ub = bam_aux2Z(ub_tag);
        const char *gene = bam_aux2Z(gx_tag);
        if (!cb || !ub || !gene) continue;

        /* seen read guard */
        uint64_t qhash = hash_qname(bam_get_qname(b), (flag & BAM_FREAD1)?1:2);
        gpointer qk = (gpointer)(uintptr_t)qhash;
        if (g_hash_table_contains(seen_reads, qk)) continue;
        g_hash_table_insert(seen_reads, qk, GUINT_TO_POINTER(1));

        int ok=0; uint64_t umi64 = pack_umi_2bit(ub,&ok);
        if (!ok) continue;

        guint32 cb_id = intern_get(cb_tab, cb);
        guint32 gene_idx = intern_get(gene_tab, gene);
        guint8 sample_idx = sample_map_get_index(smap, cb);
        uint64_t cbk = ((uint64_t)sample_idx<<56) | cb_id;

        TripKey tk = { .cbk=cbk, .gene_idx=gene_idx, .umi64=umi64 };
        TripKey *tk_heap = g_slice_new(TripKey);
        *tk_heap = tk;
        gboolean inserted = g_hash_table_insert(shd->dedup_set, tk_heap, GUINT_TO_POINTER(1));
        if (inserted) {
            uint64_t *pair = g_new(uint64_t,1);
            *pair = pair64(cb_id, gene_idx);
            gpointer val = g_hash_table_lookup(shd->counts_map, pair);
            if (!val) {
                g_hash_table_insert(shd->counts_map, pair, GUINT_TO_POINTER(1));
            } else {
                guint cnt = GPOINTER_TO_UINT(val)+1;
                g_hash_table_replace(shd->counts_map, pair, GUINT_TO_POINTER(cnt));
                g_free(pair); /* replaced key freed */
            }
        } else {
            g_slice_free(TripKey, tk_heap);
        }
        used++;
    }
    fprintf(stderr,"records=%lld processed=%lld unique_triplets=%u NNZ=%u CBs=%u Genes=%u\n",
            recs, used, g_hash_table_size(shd->dedup_set), g_hash_table_size(shd->counts_map),
            cb_tab->id_to_str->len, gene_tab->id_to_str->len);

    /* Write outputs */
    if (mkdir_p(cfg->outdir)) { fprintf(stderr,"Failed to create outdir\n"); exit(EXIT_FAILURE);}    
    char path[FILENAME_LENGTH];
    snprintf(path,sizeof(path), "%s/barcodes.tsv", cfg->outdir);
    FILE *fbc = fopen(path,"w");
    for (guint i=0;i<cb_tab->id_to_str->len;i++) fprintf(fbc, "%s\n", id_to_str(cb_tab,i));
    fclose(fbc);
    snprintf(path,sizeof(path), "%s/features.tsv", cfg->outdir);
    FILE *ffe = fopen(path,"w");
    for (guint i=0;i<gene_tab->id_to_str->len;i++) fprintf(ffe, "%s\n", id_to_str(gene_tab,i));
    fclose(ffe);
    snprintf(path,sizeof(path), "%s/matrix.mtx", cfg->outdir);
    FILE *fmtx = fopen(path,"w");
    fprintf(fmtx,"%%MatrixMarket matrix coordinate integer general\n");
    fprintf(fmtx,"%u %u %u\n", gene_tab->id_to_str->len, cb_tab->id_to_str->len, g_hash_table_size(shd->counts_map));
    /* iterate counts_map */
    GHashTableIter it; gpointer key, value;
    g_hash_table_iter_init(&it, shd->counts_map);
    while (g_hash_table_iter_next(&it, &key, &value)) {
        uint64_t pair = *(uint64_t*)key;
        guint32 cb_id = pair>>32;
        guint32 gene_idx = pair & 0xffffffffU;
        guint count = GPOINTER_TO_UINT(value);
        fprintf(fmtx, "%u %u %u\n", gene_idx+1, cb_id+1, count);
    }
    fclose(fmtx);
    write_barcode_map(cfg->outdir, cb_tab, smap);

    /* Clean */
    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    sam_close(in);
}

// MULTI-THREAD SHARD IMPLEMENTATION START

typedef struct {
    uint32_t cb_id;
    uint32_t gene_idx;
    uint64_t umi64;
    uint8_t sample_idx;
} Rec;

typedef struct {
    GAsyncQueue *queue; /* Rec* */
    shard_data  *shard;
    volatile int done;
} worker_ctx;

static gpointer shard_worker(gpointer data) {
    worker_ctx *ctx = data;
    while (1) {
        Rec *r = g_async_queue_pop(ctx->queue);
        if (!r) break; /* sentinel */
        uint64_t cbk = ((uint64_t)r->sample_idx<<56) | r->cb_id;
        TripKey tk = { .cbk = cbk, .gene_idx = r->gene_idx, .umi64 = r->umi64 };
        TripKey *heap = g_slice_new(TripKey);
        *heap = tk;
        gboolean inserted = g_hash_table_insert(ctx->shard->dedup_set, heap, GUINT_TO_POINTER(1));
        if (inserted) {
            uint64_t *pair = g_new(uint64_t,1);
            *pair = pair64(r->cb_id, r->gene_idx);
            gpointer val = g_hash_table_lookup(ctx->shard->counts_map, pair);
            if (!val) {
                g_hash_table_insert(ctx->shard->counts_map, pair, GUINT_TO_POINTER(1));
            } else {
                guint cnt = GPOINTER_TO_UINT(val)+1;
                g_hash_table_replace(ctx->shard->counts_map, pair, GUINT_TO_POINTER(cnt));
                g_free(pair);
            }
        } else {
            g_slice_free(TripKey, heap);
        }
        g_free(r);
    }
    ctx->done = 1;
    return NULL;
}

static void merge_counts(GHashTable *dest, GHashTable *src) {
    GHashTableIter it; gpointer key, val;
    g_hash_table_iter_init(&it, src);
    while (g_hash_table_iter_next(&it, &key, &val)) {
        uint64_t pair_val = *(uint64_t*)key;
        uint64_t *dup_key = g_new(uint64_t,1);
        *dup_key = pair_val;
        guint cnt_src = GPOINTER_TO_UINT(val);
        gpointer existing = g_hash_table_lookup(dest, dup_key);
        if (!existing) {
            g_hash_table_insert(dest, dup_key, GUINT_TO_POINTER(cnt_src));
        } else {
            guint cnt_total = GPOINTER_TO_UINT(existing) + cnt_src;
            g_hash_table_replace(dest, dup_key, GUINT_TO_POINTER(cnt_total));
        }
    }
}

static void process_bam_multi(cfg_t *cfg, intern_table *cb_tab, intern_table *gene_tab, sample_map *smap) {
    int N = cfg->consumer_threads;
    worker_ctx *workers = g_new0(worker_ctx, N);
    GThread **threads = g_new0(GThread*, N);
    for (int i=0;i<N;i++) {
        workers[i].queue = g_async_queue_new();
        workers[i].shard = shard_data_new();
        workers[i].done = 0;
        threads[i] = g_thread_new(NULL, shard_worker, &workers[i]);
    }

    GHashTable *seen_reads = g_hash_table_new(g_direct_hash, g_direct_equal);

    samFile *in = sam_open(cfg->bam_path, "r");
    if (!in) { perror("sam_open"); exit(EXIT_FAILURE);}    
    bam_hdr_t *hdr = sam_hdr_read(in);
    hts_set_threads(in, cfg->hts_threads);
    bam1_t *b = bam_init1();

    long long recs=0, used=0;
    while (sam_read1(in, hdr, b) >= 0) {
        if (cfg->max_records && recs>=cfg->max_records) break;
        recs++;
        uint16_t flag = b->core.flag;
        if (cfg->primary_only && (flag & (BAM_FSECONDARY|BAM_FSUPPLEMENTARY))) continue;
        if (cfg->skip_dup && (flag & BAM_FDUP)) continue;
        if (b->core.qual < cfg->min_mapq) continue;
        const uint8_t *cb_tag = bam_aux_get(b, cfg->cb_tag);
        const uint8_t *ub_tag = bam_aux_get(b, cfg->ub_tag);
        const uint8_t *gx_tag = bam_aux_get(b, cfg->gene_tag);
        if (!gx_tag) gx_tag = bam_aux_get(b, FALLBACK_GENE_TAG);
        if (!cb_tag || !ub_tag || !gx_tag) continue;
        const char *cb = bam_aux2Z(cb_tag);
        const char *ub = bam_aux2Z(ub_tag);
        const char *gene = bam_aux2Z(gx_tag);
        if (!cb || !ub || !gene) continue;
        uint64_t qhash = hash_qname(bam_get_qname(b), (flag & BAM_FREAD1)?1:2);
        gpointer qk = (gpointer)(uintptr_t)qhash;
        if (g_hash_table_contains(seen_reads, qk)) continue;
        g_hash_table_insert(seen_reads, qk, GUINT_TO_POINTER(1));
        int ok=0; uint64_t umi64 = pack_umi_2bit(ub,&ok);
        if (!ok) continue;
        guint32 cb_id = intern_get(cb_tab, cb);
        guint32 gene_idx = intern_get(gene_tab, gene);
        guint8 sample_idx = sample_map_get_index(smap, cb);
        int shard = cb_id % N;
        Rec *r = g_new(Rec,1);
        r->cb_id=cb_id; r->gene_idx=gene_idx; r->umi64=umi64; r->sample_idx=sample_idx;
        g_async_queue_push(workers[shard].queue, r);
        used++;
    }

    /* send sentinel */
    for (int i=0;i<N;i++) g_async_queue_push(workers[i].queue, NULL);
    for (int i=0;i<N;i++) g_thread_join(threads[i]);

    /* Merge shards */
    shard_data *merged = shard_data_new();
    for (int i=0;i<N;i++) {
        merge_counts(merged->counts_map, workers[i].shard->counts_map);
        // no need to merge dedup_set for final output
    }

    fprintf(stderr,"records=%lld processed=%lld NNZ=%u CBs=%u Genes=%u\n", recs, used, g_hash_table_size(merged->counts_map), cb_tab->id_to_str->len, gene_tab->id_to_str->len);

    /* outputs */
    if (mkdir_p(cfg->outdir)) { fprintf(stderr,"Failed to create outdir\n"); exit(EXIT_FAILURE);}    
    char path[FILENAME_LENGTH];
    snprintf(path,sizeof(path), "%s/barcodes.tsv", cfg->outdir);
    FILE *fbc = fopen(path,"w");
    for (guint i=0;i<cb_tab->id_to_str->len;i++) fprintf(fbc, "%s\n", id_to_str(cb_tab,i));
    fclose(fbc);
    snprintf(path,sizeof(path), "%s/features.tsv", cfg->outdir);
    FILE *ffe = fopen(path,"w");
    for (guint i=0;i<gene_tab->id_to_str->len;i++) fprintf(ffe, "%s\n", id_to_str(gene_tab,i));
    fclose(ffe);
    snprintf(path,sizeof(path), "%s/matrix.mtx", cfg->outdir);
    FILE *fmtx = fopen(path,"w");
    fprintf(fmtx,"%%MatrixMarket matrix coordinate integer general\n");
    fprintf(fmtx,"%u %u %u\n", gene_tab->id_to_str->len, cb_tab->id_to_str->len, g_hash_table_size(merged->counts_map));
    GHashTableIter it; gpointer key, val;
    g_hash_table_iter_init(&it, merged->counts_map);
    while (g_hash_table_iter_next(&it, &key, &val)) {
        uint64_t pair = *(uint64_t*)key;
        guint32 cb_id = pair>>32; guint32 gene_idx = pair & 0xffffffffU;
        guint cnt = GPOINTER_TO_UINT(val);
        fprintf(fmtx, "%u %u %u\n", gene_idx+1, cb_id+1, cnt);
    }
    fclose(fmtx);
    write_barcode_map(cfg->outdir, cb_tab, smap);

    /* TODO: free resources */
}

// MULTI-THREAD SHARD IMPLEMENTATION END

/* ---------------- barcode→sample map writer ---------------- */
static void write_barcode_map(const char *outdir, intern_table *cb_tab, sample_map *smap){
    char path[FILENAME_LENGTH];
    snprintf(path,sizeof(path), "%s/barcode_to_sample.tsv", outdir);
    FILE *f = fopen(path,"w");
    if(!f){ perror("barcode_to_sample.tsv"); return; }
    fprintf(f,"col_idx\tCB\tsample_id\n");
    /* build index→name lookup (0..255) */
    const char *names[256];
    for(int i=0;i<256;i++) names[i] = "undetermined";
    gboolean have_map = (smap && smap->sample_to_index && g_hash_table_size(smap->sample_to_index) > 0);
    if (have_map){
        GHashTableIter it; gpointer key,val;
        g_hash_table_iter_init(&it, smap->sample_to_index);
        while(g_hash_table_iter_next(&it,&key,&val)){
            guint8 idx = GPOINTER_TO_UINT(val);
            if(idx) names[idx] = (const char*)key; /* sample_id string */
        }
    }
    unsigned int n_assigned=0, n_unspecified=0, n_undetermined=0;
    for (guint i=0;i<cb_tab->id_to_str->len;i++){
        const char *cb = id_to_str(cb_tab,i);
        if (have_map){
            guint8 idx = sample_map_get_index(smap, cb);
            if (idx){
                fprintf(f, "%u\t%s\t%s\n", i+1, cb, names[idx]);
                n_assigned++;
            } else {
                fprintf(f, "%u\t%s\tunspecified\n", i+1, cb);
                n_unspecified++;
            }
        } else {
            fprintf(f, "%u\t%s\tundetermined\n", i+1, cb);
            n_undetermined++;
        }
    }
    fclose(f);
    fprintf(stderr, "CB mapping summary: assigned=%u unspecified=%u undetermined=%u\n",
            n_assigned, n_unspecified, n_undetermined);
}

int main(int argc, char **argv) {
    cfg_t cfg = {0};
    cfg.primary_only = 1; cfg.skip_dup = 1; cfg.min_mapq = 0; cfg.hts_threads=2; cfg.consumer_threads=1;
    parse_cli(&cfg, argc, argv);

    intern_table *cb_tab = intern_table_new();
    intern_table *gene_tab = intern_table_new();
    sample_map *smap = sample_map_new();
    if (cfg.sample_map_path) load_cell_sample_map(smap, cfg.sample_map_path);

    if (cfg.consumer_threads==1) {
        process_bam_single(&cfg, cb_tab, gene_tab, smap);
        write_barcode_map(cfg.outdir, cb_tab, smap);
    } else {
        process_bam_multi(&cfg, cb_tab, gene_tab, smap);
        write_barcode_map(cfg.outdir, cb_tab, smap);
    }

    return 0;
}
