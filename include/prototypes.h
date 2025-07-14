#ifndef PROTOTYPES_H
#define PROTOTYPES_H

#include "common.h"

// Function prototypes from assignBarcodes.c
void initialize_statistics(statistics *stats);
void initseq2Code();
void initcode2seq();
void initdiff2hamming(unsigned char *difference);
void initialize_complement();
feature_arrays* read_features_file(const char* filename);
unsigned char* read_whiteList(char *whitelist_filename,GHashTable *hash, int reverse_complement_flag);
void initialize_unit_sizes();
void organize_fastq_files_by_directory(int positional_arg_count, int argc, char *argv[], int optind, char *barcodeFastqFilesString, char *forwardFastqFilesString, char *reverseFastqFilesString, fastq_files_collection *fastq_files, char *barcode_pattern, char *forward_pattern, char *reverse_pattern);
void organize_fastq_files_by_type(int positional_arg_count, int argc, char *argv[], int optind, char *barcodeFastqFilesString, char *forwardFastqFilesString, char *reverseFastqFilesString, fastq_files_collection *fastq_files, char *barcode_pattern, char *forward_pattern, char *reverse_pattern, int sample_flag);
int calculate_initial_threads(fastq_files_collection *fastq_files, int parallel_by_file, int available_threads, int *consumer_threads_per_set, int *search_threads_per_consumer, int *max_concurrent_processes, int set_consumer_threads_per_set, int set_search_threads_per_consumer);
void populate_sample_args(sample_args *args, int sample_index,char *directory, fastq_files_collection *fastq_files, feature_arrays *features, int maxHammingDistance, int nThreads, memory_pool_collection *pools, statistics *stats, data_structures *hashes, uint16_t stringency, uint16_t min_counts, int barcode_constant_offset, int feature_constant_offset, int read_buffer_lines, int average_read_length, int parallel_by_file, double min_posterior, int consumer_threads_per_set);
void process_files_in_sample(sample_args *args);
void cleanup_sample(memory_pool_collection *pools, data_structures *hashes);
void free_feature_arrays(feature_arrays *features);
void initialize_data_structures(data_structures *hashes);
memory_pool_collection* initialize_memory_pool_collection();
int is_directory(const char *path);
int existing_output_skip(char keep_existing, char *directory);
double get_time_in_seconds();
int string2code(char *string, int sequence_length, unsigned char *code);
char check_sequence(char *sequence, int sequence_length);
void reverse_complement_sequence(char *sequence,  char *reverse, int length);
int generate_sequences(char* string, int string_length, int *indices, char **output, int index, int number_of_indices, char *output_indices[], int number_of_outputs);
void update_feature_counts_from_code(unsigned char *code, char *umi,  uint32_t feature_index, data_structures *hashes, memory_pool_collection *pools);
void update_umi_counts(unsigned char *code, char *umi,  uint32_t feature_index,data_structures *hashes, memory_pool_collection *pools);
void find_connected_component(gpointer start_key, uint32_t *counts, data_structures *hashes);
void add_deduped_count(feature_counts *s, uint32_t *counts, uint16_t stringency, uint16_t min_counts);
int find_variant_match(unsigned char *code, int sequence_index, unsigned char *corrected_bases );
unsigned char* find_best_posterior_match (unmatched_barcodes_features_block *entry_block, int number_of_features, double min_posterior, statistics *stats, data_structures *hashes);
void reverse_complement_in_place(char *seq);
void reverse_in_place(char *str);
void process_feature_sequence(char *sequence, feature_arrays *features, int maxHammingDistance, int nThreads, int feature_constant_offset, int max_feature_n, uint32_t *feature_index, int *hamming_distance, char *matching_sequence, uint16_t *match_position);

#endif // PROTOTYPES_H