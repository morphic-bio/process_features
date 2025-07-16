#ifndef COMMON_H
#define COMMON_H

#include <ctype.h>
#include <getopt.h>
#include <glib.h>
#include <pthread.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <zlib.h>
#include "sys/wait.h"
#include <dirent.h>
#include <stdatomic.h>
#include <sys/mman.h>
#include <errno.h>

// defines that affect structs
#define BARCODE_LENGTH 16
#define UMI_LENGTH 12
#define BARCODE_CODE_LENGTH ((BARCODE_LENGTH + 3) / 4)
#define UMI_CODE_LENGTH ((UMI_LENGTH + 3) / 4)
#define MAX_FEATURES 128
#define MAX_BARCODE_MISMATCHES 3

// other defines
#define MIN_POSTERIOR 0.975
#define MAX_FEATURE_N 3
#define MAX_BARCODE_N 1
#define MAX_MALLOCS 1000
#define MAX_FEATURE_SEQUENCE_LENGTH 64
#define READ_BUFFER_LINES 1024
#define AVERAGE_READ_LENGTH 300
#define FILENAME_LENGTH 2048
#define LINE_LENGTH 1024
#define NAME_LENGTH 40
#define MAX_FEATURE_CODE_LENGTH 40
#define FEATURE_NAME_LENGTH 32
#define LOOKUP_STRING "AAAAAAACAAAGAAATAACAAACCAACGAACTAAGAAAGCAAGGAAGTAATAAATCAATGAATTACAAACACACAGACATACCAACCCACCGACCTACGAACGCACGGACGTACTAACTCACTGACTTAGAAAGACAGAGAGATAGCAAGCCAGCGAGCTAGGAAGGCAGGGAGGTAGTAAGTCAGTGAGTTATAAATACATAGATATATCAATCCATCGATCTATGAATGCATGGATGTATTAATTCATTGATTTCAAACAACCAAGCAATCACACACCCACGCACTCAGACAGCCAGGCAGTCATACATCCATGCATTCCAACCACCCAGCCATCCCACCCCCCCGCCCTCCGACCGCCCGGCCGTCCTACCTCCCTGCCTTCGAACGACCGAGCGATCGCACGCCCGCGCGCTCGGACGGCCGGGCGGTCGTACGTCCGTGCGTTCTAACTACCTAGCTATCTCACTCCCTCGCTCTCTGACTGCCTGGCTGTCTTACTTCCTTGCTTTGAAAGAACGAAGGAATGACAGACCGACGGACTGAGAGAGCGAGGGAGTGATAGATCGATGGATTGCAAGCACGCAGGCATGCCAGCCCGCCGGCCTGCGAGCGCGCGGGCGTGCTAGCTCGCTGGCTTGGAAGGACGGAGGGATGGCAGGCCGGCGGGCTGGGAGGGCGGGGGGGTGGTAGGTCGGTGGGTTGTAAGTACGTAGGTATGTCAGTCCGTCGGTCTGTGAGTGCGTGGGTGTGTTAGTTCGTTGGTTTTAAATAACTAAGTAATTACATACCTACGTACTTAGATAGCTAGGTAGTTATATATCTATGTATTTCAATCACTCAGTCATTCCATCCCTCCGTCCTTCGATCGCTCGGTCGTTCTATCTCTCTGTCTTTGAATGACTGAGTGATTGCATGCCTGCGTGCTTGGATGGCTGGGTGGTTGTATGTCTGTGTGTTTTAATTACTTAGTTATTTCATTCCTTCGTTCTTTGATTGCTTGGTTGTTTTATTTCTTTGTTTT"
#define VISITED 255
#define CELL_SIZE 10
#define BAR_WIDTH 20
#define BASE_PADDING 10
#define BAR_GRAPH_HEIGHT 100

// Debug print function
#define DEBUG_PRINT(fmt, ...) \
    do { \
        const char *env_debug = getenv("DEBUG"); \
        if ((env_debug && atoi(env_debug) != 0) || debug) { \
            fprintf(stderr, fmt, ##__VA_ARGS__); \
        } \
    } while (0)

#define INITIAL_CAPACITY 16

// Structs
// Structs for memory management
typedef struct _storage_block {
    struct _storage_block *next;
    unsigned char *storage;
} storage_block;
typedef struct feature_arrays {
    int number_of_features;
    int max_length;
    char **feature_names;
    char *feature_names_storage;
    unsigned int *feature_lengths;
    unsigned char *feature_code_lengths;
    char **feature_sequences;
    char *feature_sequences_storage;
    unsigned char **feature_codes;
    unsigned char *feature_codes_storage;
} feature_arrays;

typedef struct feature_counts {
    unsigned char sequence_code[4];
    GHashTable *counts;
} feature_counts;

typedef struct feature_sequences {
    uint32_t counts;
    char hamming_distance;
    uint32_t feature_index;
    uint16_t match_position;
    char sequence[];
} feature_sequences;

typedef struct feature_umi_counts {
    unsigned char sequence_umi_code[8];
    GHashTable *counts;
} feature_umi_counts;

typedef struct unmatched_barcodes_features_block {
    struct unmatched_barcodes_features_block *next;
    unsigned char storage[];
} unmatched_barcodes_features_block;

typedef struct unmatched_barcodes_features {
    struct unmatched_barcodes_features_block *next;
    uint32_t feature_index;
    unsigned char *barcode;
    unsigned char *umi;
    unsigned char number_of_closest_barcodes;
    unsigned char *closest_barcodes;
    unsigned char *Qscores;
    uint16_t match_position;
} unmatched_barcodes_features;

typedef struct unmatched_barcodes_features_block_list {
    unmatched_barcodes_features_block *first_entry;
    unmatched_barcodes_features_block *last_entry;
} unmatched_barcodes_features_block_list;

typedef struct unit_sizes {
    size_t feature_counts;
    size_t feature_umi_counts;
    size_t feature_sequences;
    size_t unmatched_barcodes_features_block;
} unit_sizes;

typedef struct statistics {
    double start_time;
    size_t nMismatches;
    size_t recovered;
    size_t pending;
    size_t valid;
    size_t pending_recovered;
    size_t total_unmatched_features;
    size_t number_of_reads;
    unmatched_barcodes_features_block_list unmatched_list;
} statistics;

typedef struct _queue {
    uint64_t *data;
    size_t front;
    size_t back;
    size_t queue_size;
    size_t capacity;
} Queue;

typedef struct data_structures {
    GHashTable *filtered_hash;
    GHashTable *sequence_umi_hash;
    GHashTable *unique_features_match;
    Queue *neighbors_queue;
} data_structures;

typedef struct fastq_files_collection {
    char **barcode_fastq;
    char **forward_fastq;
    char **reverse_fastq;
    char *concatenated_files;
    int nbarcode_files;
    int nforward_files;
    int nreverse_files;
    char **sample_names;
    char *concatenated_sample_names;
    int *sample_sizes;
    int *sample_offsets;
    int nsamples;
    int max_sample_size;
    int *sorted_index;
} fastq_files_collection;

//need to have this here for sample_args
typedef struct _memory_pool {
    size_t block_size;
    size_t blocks_per_pool;
    size_t free_blocks;
    void *next_free;
    storage_block *first_block;
    storage_block *current_block;
} memory_pool;

typedef struct _memory_pool_collection {
    memory_pool *feature_counts_pool;
    memory_pool *feature_umi_counts_pool;
    memory_pool *feature_sequences_pool;
    memory_pool *unmatched_barcodes_features_block_pool;
} memory_pool_collection;

typedef struct sample_args {
    int sample_index;
    char *directory;
    fastq_files_collection *fastq_files;
    feature_arrays *features;
    int maxHammingDistance;
    int nThreads;
    memory_pool_collection *pools;
    statistics *stats;
    data_structures *hashes;
    uint16_t stringency;
    uint16_t min_counts;
    int read_buffer_lines;
    int average_read_length;
    int barcode_constant_offset;
    int feature_constant_offset;
    int parallel_by_file;
    double min_posterior;
    int consumer_threads_per_set;
} sample_args;

typedef struct fastq_reader {
    char *concatenated_filenames;
    char **filenames;
    int nfiles;
    gzFile gz_pointer;
    int filetype;
    pthread_mutex_t mutex;
    pthread_cond_t can_produce;
    pthread_cond_t can_consume;
    char **buffer;
    char *buffer_storage;
    size_t read_buffer_lines;
    size_t produce_index;
    size_t consume_index;
    size_t filled;
    int done;
} fastq_reader;

typedef struct fastq_reader_set {
    struct fastq_reader *barcode_reader;
    struct fastq_reader *forward_reader;
    struct fastq_reader *reverse_reader;
} fastq_reader_set;

typedef struct fastq_processor {
    int thread_id; // Add thread_id to identify the consumer thread
    int nsets;
    int nreaders;
    struct fastq_reader_set **reader_sets;
    struct sample_args *sample_args;
    pthread_mutex_t process_mutex;
    pthread_cond_t can_process;
} fastq_processor;

typedef struct fastq_readers_args {
    int thread_id;
    int set_index;
    int reader_type;
    int nsets;
    int nfiles;
    fastq_reader_set **reader_sets;
} fastq_readers_args;

// Queue function prototypes
void init_queue(Queue *queue);
void clear_queue(Queue *queue);
int is_empty(Queue *queue);
size_t size(Queue *queue);
void enqueue(Queue *queue, uint64_t value);
uint64_t dequeue(Queue *queue);
uint64_t peek(Queue *queue);
void free_queue(Queue *queue);

#endif // COMMON_H