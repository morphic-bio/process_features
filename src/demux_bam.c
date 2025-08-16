#include "../include/common.h"
#include "../include/utils.h"
#include "../include/globals.h"
#include "../include/io.h"
#include "../include/barcode_match.h"
#include "../include/prototypes.h"   /* for code2string() */
#include "../include/memory.h"
#include <htslib/sam.h>

static memory_pool_collection *global_pools = NULL;

#if GLIB_CHECK_VERSION(2,68,0)
#define MEMDUP g_memdup2
#else
#define MEMDUP g_memdup
#endif

#define PROBE_LEN 8

/*
 * demux_bam.c – Read STARsolo BAM, deduplicate (CB,UB,gene) and output MatrixMarket.
 * Features = sample barcode indices (probe variants). Dedup ignores sample index.
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

/* ---------------- Sample probe features ---------------- */
static feature_arrays *probe_fa = NULL;
static int direct_probe = -1; /* -1 auto; 0 hash; 1 direct */
static int probe_offset = 68;  /* default */

static inline char bam_code_to_base(uint8_t c){
    switch(c){
        case 1: return 'A';
        case 2: return 'C';
        case 4: return 'G';
        case 8: return 'T';
        default: return 'N';
    }
}

static inline char complement_base(char b)
{
    switch (b) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default : return 'N';
    }
}

static gboolean extract_kmer_from_bam(const bam1_t *b, int offset,
                                      char out8[PROBE_LEN + 1])
{
    int n = b->core.l_qseq;
    if (offset < 0 || offset + PROBE_LEN > n) return FALSE;

    const uint8_t *seq = bam_get_seq(b);

    if (b->core.flag & BAM_FREVERSE) {              /* reverse-strand read */
        int start = n - offset - PROBE_LEN;         /* 5′ on reference   */
        if (start < 0) return FALSE;
        for (int i = 0; i < PROBE_LEN; ++i) {
            /* walk right→left through the read, left→right into out8 */
            char base = bam_code_to_base(bam_seqi(seq,
                                                  start + PROBE_LEN - 1 - i));
            if (base != 'A' && base != 'C' && base != 'G' && base != 'T')
                return FALSE;
            out8[i] = complement_base(base);
        }
    } else {                                        /* forward-strand read */
        for (int i = 0; i < PROBE_LEN; ++i) {
            char base = bam_code_to_base(bam_seqi(seq, offset + i));
            if (base != 'A' && base != 'C' && base != 'G' && base != 'T')
                return FALSE;
            out8[i] = base;
        }
    }

    out8[PROBE_LEN] = '\0';
    return TRUE;
}


/* TripKey dedup helpers removed */

/* ---------------- Shard data ---------------- */

typedef struct {
    GHashTable *cb_counts;     /* key: cb_id (direct)  val: uint16_t[nprobes] */
    guint64     total_reads;
    guint64     usable_reads;
    guint64    *probe_tot;     /* per-probe total usable reads */
} shard_data;

static shard_data* shard_data_new(int nprobes){
    shard_data *s = g_new0(shard_data,1);
    s->cb_counts = g_hash_table_new(g_direct_hash, g_direct_equal);
    s->probe_tot = g_new0(guint64, nprobes);
    return s;
}


/* ---- helper: qsort comparator for barcode strings ----------------- */
static int cmp_str (const void *a, const void *b)
{
    return strcmp (*(char* const*)a, *(char* const*)b);
}
/* ------------------------------------------------------------------- */

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
    const char *sample_probes; /* path */
    int probe_offset;
    int direct_probe;
    int save_read_to_cb;
    int search_nearby; /* if non-zero, enable fallback search for probes */
} cfg_t;

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s --bam FILE --outdir DIR [options]\n"
        "Options:\n"
        "  --cell_sample_map FILE   Optional CB→sample mapping TSV (unused in probe mode)\n"
        "  --sample_probes FILE     Probe 8-mer table (variant,unused,BCname)\n"
        "  --probe_offset N         0-based probe offset (default 68)\n"
        "  --cb_tag STR             Default CB\n"
        "  --ub_tag STR             Default UB\n"
        "  --gene_tag STR           Default GX (fallback GE)\n"
        "  --hts_threads N          BGZF threads (default 2) or -S N\n"
        "  --threads N              Consumer shards (default 1) or -t N\n"
        "  --min_mapq N             Minimum MAPQ (default 0)\n"
        "  --no_primary_filter      Keep secondary/supplementary\n"
        "  --keep_dup               Keep duplicate reads\n"
        "  --max_records N          Process first N records only\n"
        "  -v                       Debug\n",
        prog);
}

static feature_arrays* load_probe_variants_to_features(const char *path) {
    /* Three-column TSV: VARIANT_8MER  CANONICAL_8MER  BARCODE_ID */
    FILE *fp = fopen(path, "r");
    if (!fp) { perror("sample_probes"); exit(EXIT_FAILURE); }

    /* first pass – collect unique barcode IDs and their canonical sequence */
    GHashTable *bc2canon = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);
    char line[LINE_LENGTH];
    int uniq_cnt = 0, name_size = 0, seq_size = 0, code_size = 0, variant_cnt = 0;
    while (fgets(line, sizeof(line), fp)) {
        char *save=NULL;
        char *variant = strtok_r(line, "\t \r\n", &save);
        char *canonical = strtok_r(NULL, "\t \r\n", &save);
        char *bc_id   = strtok_r(NULL, "\t \r\n", &save);
        if (!variant || !bc_id) continue;
        if ((int)strlen(variant) != PROBE_LEN) continue;
        if (!check_sequence(variant, PROBE_LEN)) continue;
        if (!canonical || strlen(canonical)!=PROBE_LEN || !check_sequence(canonical, PROBE_LEN))
            canonical = variant; /* fallback */
        if (!g_hash_table_contains(bc2canon, bc_id)) {
            g_hash_table_insert(bc2canon, g_strdup(bc_id), g_strdup(canonical));
            uniq_cnt++;
            name_size += (int)strlen(bc_id) + 1;
            seq_size  += PROBE_LEN + 1;
            code_size += (PROBE_LEN+3)/4;
        }
        if (canonical != variant) variant_cnt++;    
    }

    feature_arrays *fa = allocate_feature_arrays(name_size, seq_size, code_size, uniq_cnt, PROBE_LEN);
    fa->common_length = PROBE_LEN;
    //allocate a array of blocks of size (PROBE_LEN+3)/4
    uint8_t *variant_codes = malloc(variant_cnt * (PROBE_LEN+3)/4);
    memset(variant_codes, 0, variant_cnt * (PROBE_LEN+3)/4);

    /* ---------- collect & sort unique BC-IDs ---------- */
    char **bc_list = g_new(char*, uniq_cnt);
    {
        GHashTableIter it;
        gpointer key, value;
        g_hash_table_iter_init (&it, bc2canon);
        guint i = 0;
        while (g_hash_table_iter_next (&it, &key, &value))
            bc_list[i++] = (char*)key;          /* key already dup’ed */
    }

    /* alphabetic sort */
    qsort (bc_list, uniq_cnt, sizeof(char*), cmp_str);

    /* second pass – build arrays */
    GHashTable *bc2idx = g_hash_table_new_full(g_str_hash, g_str_equal, NULL, NULL);
    int idx = 0;
    fa->feature_names[0]     = fa->feature_names_storage;
    fa->feature_sequences[0] = fa->feature_sequences_storage;
    fa->feature_codes[0]     = fa->feature_codes_storage;
    for (idx = 0; idx < uniq_cnt; idx++) {
        const char *bc_id = bc_list[idx];
        const char *canon = g_hash_table_lookup (bc2canon, bc_id);
        /* store name */
        strcpy(fa->feature_names[idx], bc_id);
        if (idx+1 < uniq_cnt)
            fa->feature_names[idx+1] = fa->feature_names[idx] + strlen(bc_id) + 1;
        /* store canonical sequence */
        strcpy(fa->feature_sequences[idx], canon);
        fa->feature_lengths[idx] = PROBE_LEN;
        /* cast away const because string2code expects char* */
        fa->feature_code_lengths[idx] =
            string2code((char *)canon, PROBE_LEN, fa->feature_codes[idx]);
        if (idx+1 < uniq_cnt) {
            fa->feature_sequences[idx+1] = fa->feature_sequences[idx] + PROBE_LEN + 1;
            fa->feature_codes[idx+1]     = fa->feature_codes[idx] + fa->feature_code_lengths[idx];
        }
        /* insert canonical code into lookup */
        if (fa->feature_lengths[idx] == PROBE_LEN) {
            GBytes *k = g_bytes_new_static(fa->feature_codes[idx], fa->feature_code_lengths[idx]);
            g_hash_table_insert(feature_code_hash, k, GUINT_TO_POINTER(idx+1));
        }
        g_hash_table_insert(bc2idx, (gpointer)bc_id, GUINT_TO_POINTER(idx + 1));
    }
    g_free(bc_list);

    /* third pass – map every variant 8-mer to the same index */
    fseek(fp, 0, SEEK_SET);
    int variant_idx = 0;
    while (fgets(line, sizeof(line), fp)) {
        char *save=NULL;
        char *variant = strtok_r(line, "\t \r\n", &save);
        char *canonical = strtok_r(NULL, "\t \r\n", &save);
        char *bc_id   = strtok_r(NULL, "\t \r\n", &save);
        if (!variant || !bc_id) continue;
        if (strcmp(variant, canonical) == 0) continue;
        if ((int)strlen(variant) != PROBE_LEN) continue;
        if (!check_sequence(variant, PROBE_LEN)) continue;

        gpointer iptr = g_hash_table_lookup(bc2idx, bc_id);
        if (!iptr) continue;              /* should not happen */
        guint bc_index = GPOINTER_TO_UINT(iptr) - 1;   /* back to 0-based */
        uint8_t *codebuf = variant_codes + variant_idx * (PROBE_LEN+3)/4;
        int codelen = string2code(variant, PROBE_LEN, codebuf);
        GBytes *k = g_bytes_new_static(codebuf, codelen);
        g_hash_table_insert(feature_code_hash, k, GUINT_TO_POINTER(bc_index+1));
        variant_idx++;
    }
    //print the size of the feature_code_hash
    fprintf(stderr, "DEBUG: variant_idx=%d\n", variant_idx);
    fprintf(stderr, "DEBUG: feature_code_hash size=%u\n", g_hash_table_size(feature_code_hash));   
    fclose(fp);
    g_hash_table_destroy(bc2canon);
    g_hash_table_destroy(bc2idx);

    fprintf(stderr, "Loaded %d unique sample barcode IDs (len=%d)\n", uniq_cnt, PROBE_LEN);
    return fa;
}

static void parse_cli(cfg_t *cfg, int argc, char **argv) {
    static struct option long_opts[] = {
        {"bam", required_argument, 0, 0},
        {"outdir", required_argument,0,0},
        {"cell_sample_map", required_argument,0,0},
        {"sample_probes", required_argument,0,0},
        {"probe_offset", required_argument,0,0},
        {"cb_tag", required_argument,0,0},
        {"ub_tag", required_argument,0,0},
        {"gene_tag", required_argument,0,0},
        {"hts_threads", required_argument,0,0},
        {"threads", required_argument,0,0},
        {"min_mapq", required_argument,0,0},
        {"no_primary_filter", no_argument,0,0},
        {"keep_dup", no_argument,0,0},
        {"max_records", required_argument,0,0},
        {"save_read_to_cb", no_argument,0,0},
        {"search_nearby", no_argument, 0, 0},
        {"debug", no_argument,0,'v'},
        {0,0,0,0}
    };
    int opt_idx=0;
    while (1) {
        int c = getopt_long(argc, argv, "vS:t:", long_opts, &opt_idx);
        if (c==-1) break;
        if (c=='v') { cfg->debug=1; continue; }
        if (c=='S') { cfg->hts_threads = atoi(optarg); continue; }
        if (c=='t') { cfg->consumer_threads = atoi(optarg); continue; }
        if (c!=0) continue;
        const char *opt_name = long_opts[opt_idx].name;
        if (!strcmp(opt_name, "bam")) cfg->bam_path = optarg;
        else if (!strcmp(opt_name, "outdir")) cfg->outdir = optarg;
        else if (!strcmp(opt_name, "cell_sample_map")) cfg->sample_map_path = optarg;
        else if (!strcmp(opt_name, "sample_probes")) cfg->sample_probes = optarg;
        else if (!strcmp(opt_name, "probe_offset")) cfg->probe_offset = atoi(optarg);
        else if (!strcmp(opt_name, "cb_tag")) cfg->cb_tag = optarg;
        else if (!strcmp(opt_name, "ub_tag")) cfg->ub_tag = optarg;
        else if (!strcmp(opt_name, "gene_tag")) cfg->gene_tag = optarg;
        else if (!strcmp(opt_name, "hts_threads")) cfg->hts_threads = atoi(optarg);
        else if (!strcmp(opt_name, "threads")) cfg->consumer_threads = atoi(optarg);
        else if (!strcmp(opt_name, "min_mapq")) cfg->min_mapq = atoi(optarg);
        else if (!strcmp(opt_name, "no_primary_filter")) cfg->primary_only = 0;
        else if (!strcmp(opt_name, "keep_dup")) cfg->skip_dup = 0;
        else if (!strcmp(opt_name, "max_records")) cfg->max_records = atoll(optarg);
        else if (!strcmp(opt_name, "save_read_to_cb")) cfg->save_read_to_cb = 1;
        else if (!strcmp(opt_name, "search_nearby")) cfg->search_nearby = 1;
        else {
            fprintf(stderr, "Error: unsupported option --%s\n", opt_name);
            usage(argv[0]);
            exit(EXIT_FAILURE);
        }
    }
    if (!cfg->bam_path || !cfg->outdir) { usage(argv[0]); exit(EXIT_FAILURE);}    
    if (!cfg->cb_tag) cfg->cb_tag = DEFAULT_CB_TAG;
    if (!cfg->ub_tag) cfg->ub_tag = DEFAULT_UB_TAG;
    if (!cfg->gene_tag) cfg->gene_tag = DEFAULT_GENE_TAG;
    if (cfg->hts_threads<=0) cfg->hts_threads=2;
    if (cfg->consumer_threads<=0) cfg->consumer_threads=1;
    if (cfg->probe_offset<=0) cfg->probe_offset = probe_offset;
}

/* ---------------- Worker function (single-thread) ---------------- */

static void process_bam_single(cfg_t *cfg,
                               intern_table *cb_tab,
                               intern_table *gene_tab)
{
    samFile *in = sam_open(cfg->bam_path, "r");
    if (!in) { perror("sam_open"); exit(EXIT_FAILURE);}    
    bam_hdr_t *hdr = sam_hdr_read(in);
    hts_set_threads(in, cfg->hts_threads);

    int nprobes = probe_fa ? probe_fa->number_of_features : 0;
    shard_data *shd = shard_data_new(nprobes);
    /* seen_reads guard removed: we rely on primary-only filter and TripKey dedup */
    /* (no separate hash table for read IDs anymore) */
    GHashTable *read_map = NULL;
    if (cfg->save_read_to_cb)
        read_map = g_hash_table_new_full(g_str_hash, g_str_equal, free, free);

    bam1_t *b = bam_init1();
    long long recs=0;
    while (sam_read1(in, hdr, b) >= 0) {
        if (cfg->max_records && recs>=cfg->max_records) break;
        recs++;
        const char *qname = bam_get_qname(b);
        uint16_t flag = b->core.flag;
        if (cfg->primary_only && (flag & (BAM_FSECONDARY|BAM_FSUPPLEMENTARY))) {
            continue; 
        }
        if (cfg->skip_dup && (flag & BAM_FDUP)) { 
            continue; 
        }
        if (b->core.qual < cfg->min_mapq) { 
            continue; 
        }

        /* Tags */
        const uint8_t *cb_tag = bam_aux_get(b, cfg->cb_tag);
        const uint8_t *ub_tag = bam_aux_get(b, cfg->ub_tag);
        const uint8_t *gx_tag = bam_aux_get(b, cfg->gene_tag);
        if (!gx_tag) gx_tag = bam_aux_get(b, FALLBACK_GENE_TAG);
        if (!cb_tag || !ub_tag || !gx_tag) { 
            continue; 
        }
        const char *cb = bam_aux2Z(cb_tag);
        const char *ub = bam_aux2Z(ub_tag);
        const char *gene = bam_aux2Z(gx_tag);
        if (!cb || !ub || !gene) { 
            continue; 
        }

        /* packed UMI stored in 'out' but no longer used for dedup */

        /* sample probe lookup */
        int sample_idx = 0;
        if (probe_fa) {
            char kmer[PROBE_LEN + 1];
            /* try primary offset first */
            if (extract_kmer_from_bam(b, cfg->probe_offset, kmer)) {
                int idx = feature_lookup_kmer(kmer, PROBE_LEN, probe_fa, direct_probe);
                if (idx > 0 && idx <= probe_fa->number_of_features)
                    sample_idx = idx;
            }
            /* fallback search if enabled and not found */
            if (sample_idx == 0 && cfg->search_nearby) {
                int offs[4] = { cfg->probe_offset + 1,
                                cfg->probe_offset + 2,
                                cfg->probe_offset - 1,
                                cfg->probe_offset - 2 };
                for (int oi = 0; oi < 4 && sample_idx == 0; ++oi) {
                    if (!extract_kmer_from_bam(b, offs[oi], kmer))
                        continue; /* invalid offset or non-ACGT base */
                    int idx = feature_lookup_kmer(kmer, PROBE_LEN, probe_fa, direct_probe);
                    if (idx > 0 && idx <= probe_fa->number_of_features)
                        sample_idx = idx;   /* early exit once found */
                }
            }
        }

        guint32 cb_id = intern_get(cb_tab, cb);

        if (sample_idx>0){
            /* get / create per-CB counts array */
            gpointer pv = g_hash_table_lookup(shd->cb_counts, GUINT_TO_POINTER(cb_id));
            uint16_t *arr;
            if (!pv){
                /* allocate from global pools if available, else g_malloc */
                #ifdef CB_COUNTS_BLOCK_SIZE
                if (global_pools && global_pools->cb_counts_pool){
                    arr = allocate_memory_from_pool(global_pools->cb_counts_pool);
                } else {
                    arr = g_malloc0(nprobes * sizeof(uint16_t));
                }
                #else
                arr = g_malloc0(nprobes * sizeof(uint16_t));
                #endif
                g_hash_table_insert(shd->cb_counts, GUINT_TO_POINTER(cb_id), arr);
            } else {
                arr = (uint16_t*)pv;
            }
            if (arr[sample_idx-1] < UINT16_MAX)
                arr[sample_idx-1]++;

            shd->usable_reads++;
            shd->probe_tot[sample_idx-1]++;
        }
        if (sample_idx>0 && cfg->save_read_to_cb){
            char *val = g_strdup_printf("%s\t%s\t%s", cb, ub, gene);
            g_hash_table_replace(read_map, g_strdup(qname), val);
        }
        shd->total_reads++;
    }

    /* TripKey dedupulation disabled (Stage 3). Stats simplified */
    fprintf(stderr, "records=%lu usable=%lu\n",
            (unsigned long)shd->total_reads,
            (unsigned long)shd->usable_reads);


    /* Write outputs */
    if (mkdir_p(cfg->outdir)) { fprintf(stderr,"Failed to create outdir\n"); exit(EXIT_FAILURE);}    
    char path[FILENAME_LENGTH];

    /* ---- Write barcodes.tsv & build CB→column map ---- */
    GHashTableIter bit; gpointer bk,bv;
    GHashTable *cb_to_col = g_hash_table_new(g_direct_hash, g_direct_equal);
    snprintf(path,sizeof(path), "%s/barcodes.tsv", cfg->outdir);
    FILE *fbc = fopen(path,"w");
    guint32 col=0;
    g_hash_table_iter_init(&bit, shd->cb_counts);
    while (g_hash_table_iter_next(&bit, &bk, &bv)){
        guint32 cb_id = GPOINTER_TO_UINT(bk);
        fprintf(fbc, "%s\n", id_to_str(cb_tab, cb_id));
        g_hash_table_insert(cb_to_col, bk, GUINT_TO_POINTER(++col));
    }
    fclose(fbc);

    /* features.tsv = probe names */
    snprintf(path,sizeof(path), "%s/features.tsv", cfg->outdir);
    FILE *ffe = fopen(path,"w");
    for (int i=0;i<nprobes;i++) fprintf(ffe, "%s\n", probe_fa->feature_names[i]);
    fclose(ffe);

    /* ---- Compute NNZ ---- */
    guint64 nnz=0;
    g_hash_table_iter_init(&bit, shd->cb_counts);
    while (g_hash_table_iter_next(&bit, &bk, &bv)){
        uint16_t *arr = (uint16_t*)bv;
        for (int p=0;p<nprobes;p++) if (arr[p]) nnz++;
    }

    /* matrix.mtx */
    snprintf(path,sizeof(path), "%s/matrix.mtx", cfg->outdir);
    FILE *fmtx = fopen(path,"w");
    fprintf(fmtx, "%%MatrixMarket matrix coordinate integer general\n");
    fprintf(fmtx, "%d %u %llu\n", nprobes, col, (unsigned long long)nnz);

    g_hash_table_iter_init(&bit, shd->cb_counts);
    while (g_hash_table_iter_next(&bit, &bk, &bv)){
        uint32_t col_idx = GPOINTER_TO_UINT(g_hash_table_lookup(cb_to_col, bk));
        uint16_t *arr = (uint16_t*)bv;
        for (int p=0;p<nprobes;p++)
            if (arr[p]) fprintf(fmtx, "%d %u %u\n", p+1, col_idx, arr[p]);
    }
    fclose(fmtx);

    /* stats.txt */
    snprintf(path,sizeof(path), "%s/stats.txt", cfg->outdir);
    FILE *fst = fopen(path,"w");
    fprintf(fst,"Total_records\t%lu\n", (unsigned long)shd->total_reads);
    fprintf(fst,"Usable_records\t%lu\n", (unsigned long)shd->usable_reads);
    fprintf(fst,"BC_ID\tTotal_reads\n");
    for(int i=0;i<nprobes;i++){
        fprintf(fst,"%s\t%lu\n", probe_fa->feature_names[i],
                (unsigned long)shd->probe_tot[i]);
    }
    fclose(fst);

    /* Clean */
    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    sam_close(in);

    /* probe_single removed */

    /* write read_to_cb and rejected_qnames */
    if (cfg->save_read_to_cb){
        char path[FILENAME_LENGTH];
        snprintf(path, sizeof(path), "%s/read_to_cb_umi_gene.txt", cfg->outdir);
        FILE *fp = fopen(path, "w");
        if (!fp) { perror("read_to_cb"); exit(EXIT_FAILURE);} 
        GHashTableIter iter; gpointer k,v; g_hash_table_iter_init(&iter, read_map);
        while (g_hash_table_iter_next(&iter, &k, &v)) fprintf(fp, "%s\t%s\n", (char*)k, (char*)v);
        fclose(fp);
        g_hash_table_destroy(read_map);
    }
}

int main(int argc, char **argv) {
    cfg_t cfg = {0};
    cfg.primary_only = 1; cfg.skip_dup = 1; cfg.min_mapq = 0; cfg.hts_threads=2; cfg.consumer_threads=1;
    cfg.probe_offset = 68; cfg.direct_probe = -1;
    parse_cli(&cfg, argc, argv);

    /* init matching tables */
    barcode_match_init();
    feature_code_hash = g_hash_table_new_full(g_bytes_hash,
                                              g_bytes_equal,
                                              (GDestroyNotify) g_bytes_unref,
                                              NULL);

    if (cfg.sample_probes) {
        probe_fa     = load_probe_variants_to_features(cfg.sample_probes);
        probe_offset = cfg.probe_offset;
        direct_probe = 0;          /* disable direct search for now */
    }

    /* ---- Stage 4: set up memory pools for CB×probe counts ---- */
    dynamic_struct_sizes.cb_probe_counts = (probe_fa?probe_fa->number_of_features:1) * sizeof(uint16_t);
    global_pools = initialize_memory_pool_collection();

    intern_table *cb_tab = intern_table_new();
    intern_table *gene_tab = intern_table_new();

    process_bam_single(&cfg, cb_tab, gene_tab);

    if (global_pools) free_memory_pool_collection(global_pools);
    return 0;
}
