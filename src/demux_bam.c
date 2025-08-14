#include "../include/common.h"
#include "../include/utils.h"
#include "../include/globals.h"
#include "../include/io.h"
#include "../include/barcode_match.h"
#include "../include/prototypes.h"   /* for code2string() */
#include <htslib/sam.h>

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

/* ------------ DEBUG: dump all feature-hash keys -------------- */
static void dump_feature_codes(void)
{
    GHashTableIter it;
    gpointer key, val;
    g_hash_table_iter_init(&it, feature_code_hash);

    char seqbuf[PROBE_LEN * 4 + 1];   /* enough for 4 bases/byte */
    while (g_hash_table_iter_next(&it, &key, &val)) {
        GBytes *g = (GBytes *)key;
        gsize blen;
        const guint8 *data = g_bytes_get_data(g, &blen);

        /* convert the packed code back to a string */
        code2string((unsigned char *)data, seqbuf, (int)blen);

        /* trim to the real probe length when printing */
        seqbuf[PROBE_LEN] = '\0';
        fprintf(stderr, "[feature_hash] %s → %u\n",
                seqbuf, GPOINTER_TO_UINT(val));
    }
    fprintf(stderr, "[feature_hash] ---- end dump ----\n");
}
/* ------------------------------------------------------------- */

/* ---------------- TripKey hashing ---------------- */

typedef struct {
    uint64_t cbk;       /* sample assignment not used; lower 56 bits hold cb_id */
    uint32_t gene_idx;  /* used for dedup only */
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

/* Key for counts (cb_id,sample_idx) packed into uint64 */
static inline uint64_t pair64(uint32_t a, uint32_t b) {
    return ((uint64_t)a<<32)|b;
}

/* ---------------- Shard data ---------------- */

typedef struct {
    GHashTable *dedup_set;     /* TripKey* -> votes (GHashTable<direct-int key -> count>) */
    GHashTable *counts_map;    /* uint64_t* -> uint32 count (GUINT_TO_POINTER) */
    GHashTable *ambiguous_cb;  /* direct-int key cb_id -> 1 */
    guint64      total_reads;
    guint64      usable_reads;
    guint64      single_reads;
    guint64      multi_reads;
    guint64     *probe_tot;    /* per-probe total usable reads */
    guint64     *probe_single; /* per-probe reads from single-barcode cells */
} shard_data;

static void tripkey_free(gpointer data){ g_slice_free(TripKey, data); }

static shard_data* shard_data_new(int nprobes) {
    shard_data *s = g_new0(shard_data,1);
    s->dedup_set = g_hash_table_new_full(tripkey_hash, tripkey_equal, tripkey_free, (GDestroyNotify)g_hash_table_destroy);
    s->counts_map = g_hash_table_new_full(g_int64_hash, g_int64_equal, g_free, NULL);
    s->ambiguous_cb = g_hash_table_new(g_direct_hash, g_direct_equal);
    s->probe_tot    = g_new0(guint64, nprobes);
    s->probe_single = g_new0(guint64, nprobes);
    return s;
}

/* ---------------- Seen reads table ---------------- */
static guint64 hash_qname(const char *qname, int read12) {
    uint64_t h=14695981039346656037ULL;
    for (const unsigned char *p=(const unsigned char*)qname; *p; ++p) { h ^= *p; h *= 1099511628211ULL; }
    h ^= read12;
    return h;
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
        fprintf(stderr, "DEBUG: %d variant=%s\tbc_id=%s\n", variant_idx, variant, bc_id);
        uint8_t *codebuf = variant_codes + variant_idx * (PROBE_LEN+3)/4;
        int codelen = string2code(variant, PROBE_LEN, codebuf);
        GBytes *k = g_bytes_new_static(codebuf, codelen);
        g_hash_table_insert(feature_code_hash, k, GUINT_TO_POINTER(bc_index+1));
        variant_idx++;
    }
    //print the size of the feature_code_hash
    fprintf(stderr, "DEBUG: variant_idx=%d\n", variant_idx);
    fprintf(stderr, "DEBUG: feature_code_hash size=%zu\n", g_hash_table_size(feature_code_hash));   
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

    /* one-time debug dump of all feature codes */
    static int dumped = 0;
    if (!dumped) {
        dump_feature_codes();
        dumped = 1;
    }

    int nprobes = probe_fa ? probe_fa->number_of_features : 0;
    shard_data *shd = shard_data_new(nprobes);
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

        /* UMI pack */
        int ok=1; uint64_t out=0; int len=0; 
        for (; ub[len] && len<30; ++len){ 
            char c=ub[len]; 
            uint8_t bits=(c=='A'||c=='a')?0:(c=='C'||c=='c')?1:(c=='G'||c=='g')?2:(c=='T'||c=='t')?3:4; 
            if (bits>3){ ok=0; break;} 
            out |= ((uint64_t)bits)<<(len*2);
        } 
        out |= ((uint64_t)len)<<60; 
        if (!ok) continue; 
        uint64_t umi64 = out;

        /* sample probe lookup */
        int sample_idx = 0;
        if (probe_fa){
            char kmer[PROBE_LEN+1];
            if (extract_kmer_from_bam(b, cfg->probe_offset, kmer)){
                int idx = feature_lookup_kmer(kmer, PROBE_LEN, probe_fa, direct_probe);
                if (idx>0 && idx<=probe_fa->number_of_features) sample_idx = idx;
            }
        }

        guint32 cb_id = intern_get(cb_tab, cb);
        guint32 gene_idx = intern_get(gene_tab, gene); /* for dedup only */
        uint64_t cbk = cb_id; /* top bits unused */

        TripKey tk = { .cbk=cbk, .gene_idx=gene_idx, .umi64=umi64 };
        gpointer orig_key=NULL, value=NULL;
        if (!g_hash_table_lookup_extended(shd->dedup_set, &tk, &orig_key, &value)){
            TripKey *heap = g_slice_new(TripKey); *heap = tk;
            GHashTable *votes = g_hash_table_new(g_direct_hash, g_direct_equal);
            g_hash_table_insert(votes, GUINT_TO_POINTER(sample_idx), GUINT_TO_POINTER(1));
            g_hash_table_insert(shd->dedup_set, heap, votes);
        } else {
            GHashTable *votes = (GHashTable*)value;
            gpointer v = g_hash_table_lookup(votes, GUINT_TO_POINTER(sample_idx));
            guint cnt = v ? GPOINTER_TO_UINT(v)+1 : 1;
            g_hash_table_replace(votes, GUINT_TO_POINTER(sample_idx), GUINT_TO_POINTER(cnt));
        }
        if (sample_idx>0){
            used++;
            shd->usable_reads++;
            shd->probe_tot[sample_idx-1]++;
        }
        shd->total_reads++;
    }

    /* Majority resolution and counts */
    GHashTableIter dit; gpointer key,value;
    g_hash_table_iter_init(&dit, shd->dedup_set);
    while (g_hash_table_iter_next(&dit, &key, &value)){
        TripKey *tk = (TripKey*)key;
        GHashTable *votes = (GHashTable*)value;
        guint best_idx=0, best_cnt=0, second_idx=0, second_cnt=0;
        GHashTableIter vit; gpointer k2,v2;
        g_hash_table_iter_init(&vit, votes);
        while (g_hash_table_iter_next(&vit, &k2, &v2)){
            guint sidx = GPOINTER_TO_UINT(k2); guint cnt = GPOINTER_TO_UINT(v2);
            if (cnt > best_cnt){ second_cnt=best_cnt; second_idx=best_idx; best_cnt=cnt; best_idx=sidx; }
            else if (cnt > second_cnt){ second_cnt=cnt; second_idx=sidx; }
        }
        if (best_cnt>0 && second_cnt>0 && second_cnt==best_cnt && best_idx!=0 && second_idx!=0){
            /* tie between non-zero sample indices → drop CB */
            guint32 cb_id = (guint32)(tk->cbk & 0xffffffffULL);
            g_hash_table_insert(shd->ambiguous_cb, GUINT_TO_POINTER(cb_id), GUINT_TO_POINTER(1));
            continue;
        }
        /* accumulate single/multi stats */
        if (best_idx>0){
            if (second_cnt==0){
                shd->single_reads++;
                shd->probe_single[best_idx-1]++;
            } else {
                shd->multi_reads++;
            }
        }

        /* credit 1 to (cb_id, best_idx) if CB not ambiguous */
        guint32 cb_id = (guint32)(tk->cbk & 0xffffffffULL);
        if (g_hash_table_contains(shd->ambiguous_cb, GUINT_TO_POINTER(cb_id))) continue;
        uint64_t *pair = g_new(uint64_t,1);
        *pair = pair64(cb_id, best_idx);
        gpointer val = g_hash_table_lookup(shd->counts_map, pair);
        if (!val)
            g_hash_table_insert(shd->counts_map, pair, GUINT_TO_POINTER(1));
        else {
            guint cnt = GPOINTER_TO_UINT(val) + 1;
            g_hash_table_replace(shd->counts_map, pair, GUINT_TO_POINTER(cnt));
        }
    }

    fprintf(stderr,"records=%lu usable=%lu single=%lu multi=%lu triplets=%u NNZ=%u CBs=%u Probes=%d\n",
            (unsigned long)shd->total_reads, (unsigned long)shd->usable_reads,
            (unsigned long)shd->single_reads, (unsigned long)shd->multi_reads,
            g_hash_table_size(shd->dedup_set), g_hash_table_size(shd->counts_map),
            cb_tab->id_to_str->len, probe_fa?probe_fa->number_of_features:0);

    /* Write outputs */
    if (mkdir_p(cfg->outdir)) { fprintf(stderr,"Failed to create outdir\n"); exit(EXIT_FAILURE);}    
    char path[FILENAME_LENGTH];

    /* barcodes.tsv excluding ambiguous */
    GHashTable *cb_to_col = g_hash_table_new(g_direct_hash, g_direct_equal);
    snprintf(path,sizeof(path), "%s/barcodes.tsv", cfg->outdir);
    FILE *fbc = fopen(path,"w");
    guint32 col=0;
    for (guint i=0;i<cb_tab->id_to_str->len;i++){
        if (g_hash_table_contains(shd->ambiguous_cb, GUINT_TO_POINTER(i))) continue;
        const char *cb = id_to_str(cb_tab,i);
        fprintf(fbc, "%s\n", cb);
        g_hash_table_insert(cb_to_col, GUINT_TO_POINTER(i), GUINT_TO_POINTER(++col));
    }
    fclose(fbc);

    /* features.tsv = probe names */
    snprintf(path,sizeof(path), "%s/features.tsv", cfg->outdir);
    FILE *ffe = fopen(path,"w");
    /* nprobes already declared earlier */
    for (int i=0;i<nprobes;i++) fprintf(ffe, "%s\n", probe_fa->feature_names[i]);
    fclose(ffe);

    /* matrix.mtx: rows=probes, cols=kept CBs */
    snprintf(path,sizeof(path), "%s/matrix.mtx", cfg->outdir);
    FILE *fmtx = fopen(path,"w");
    fprintf(fmtx,"%%MatrixMarket matrix coordinate integer general\n");
    fprintf(fmtx,"%d %u %u\n", nprobes, col, g_hash_table_size(shd->counts_map));
    GHashTableIter it; gpointer k3, v3;
    g_hash_table_iter_init(&it, shd->counts_map);
    while (g_hash_table_iter_next(&it, &k3, &v3)){
        uint64_t pair = *(uint64_t*)k3;
        guint32 cb_id = pair>>32; guint32 sidx = pair & 0xffffffffU;
        gpointer colp = g_hash_table_lookup(cb_to_col, GUINT_TO_POINTER(cb_id));
        if (!colp) continue; /* ambiguous dropped */
        fprintf(fmtx, "%u %u %u\n", sidx, GPOINTER_TO_UINT(colp), GPOINTER_TO_UINT(v3));
    }
    fclose(fmtx);

    /* stats.txt */
    snprintf(path,sizeof(path), "%s/stats.txt", cfg->outdir);
    FILE *fst = fopen(path,"w");
    fprintf(fst,"Total_records\t%lu\n", (unsigned long)shd->total_reads);
    fprintf(fst,"Usable_records\t%lu\n", (unsigned long)shd->usable_reads);
    fprintf(fst,"Fraction_single_barcode\t%.6f\n", shd->usable_reads? (double)shd->single_reads/shd->usable_reads:0.0);
    fprintf(fst,"Fraction_multi_barcode\t%.6f\n", shd->usable_reads? (double)shd->multi_reads/shd->usable_reads:0.0);
    fprintf(fst,"BC_ID\tTotal_reads\tSingle_barcode_reads\n");
    for(int i=0;i<nprobes;i++){
        fprintf(fst,"%s\t%lu\t%lu\n", probe_fa->feature_names[i],
                (unsigned long)shd->probe_tot[i], (unsigned long)shd->probe_single[i]);
    }
    fclose(fst);

    /* Clean */
    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    sam_close(in);

    g_free(shd->probe_tot);
    g_free(shd->probe_single);
}

int main(int argc, char **argv) {
    cfg_t cfg = {0};
    cfg.primary_only = 1; cfg.skip_dup = 1; cfg.min_mapq = 0; cfg.hts_threads=2; cfg.consumer_threads=1;
    cfg.probe_offset = 68; cfg.direct_probe = -1;
    parse_cli(&cfg, argc, argv);

    /* init matching tables */
    barcode_match_init();
    feature_code_hash = g_hash_table_new_full(g_bytes_hash, g_bytes_equal, (GDestroyNotify)g_bytes_unref, NULL);
    if (cfg.sample_probes){ probe_fa = load_probe_variants_to_features(cfg.sample_probes); direct_probe = cfg.direct_probe; probe_offset = cfg.probe_offset; }

    intern_table *cb_tab = intern_table_new();
    intern_table *gene_tab = intern_table_new();

    process_bam_single(&cfg, cb_tab, gene_tab);

    return 0;
}
