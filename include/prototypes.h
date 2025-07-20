#ifndef PROTOTYPES_H
#define PROTOTYPES_H

#include "../include/common.h"

// Function prototypes from assignBarcodes.c
void add_deduped_count(GHashTable* temp_deduped_hash, uint32_t *clique_counts, uint16_t stringency, uint16_t min_counts);
int calculate_initial_threads(fastq_files_collection *fastq_files, int available_threads, int *consumer_threads_per_set, int *search_threads_per_consumer, int *max_concurrent_processes, int set_consumer_threads_per_set, int set_search_threads_per_consumer);
char check_sequence(char *sequence, int sequence_length);
void cleanup_sample(memory_pool_collection *pools, data_structures *hashes);
int existing_output_skip(char keep_existing, char *directory);
void find_connected_component(gpointer start_key, uint32_t *counts, data_structures *hashes);
void find_deduped_counts(data_structures *hashes, GHashTable* barcode_to_deduped_counts, uint16_t stringency, uint16_t min_counts);
int find_variant_match(unsigned char *code, int sequence_index, unsigned char *corrected_bases );
unsigned char* find_best_posterior_match (unmatched_barcodes_features_block *entry_block, int number_of_features, double min_posterior, statistics *stats, data_structures *hashes);
void free_feature_arrays(feature_arrays *features);
int generate_sequences(char* string, int string_length, int *indices, char **output, int index, int number_of_indices, char *output_indices[], int number_of_outputs);
double get_time_in_seconds();
void initcode2seq();
void initialize_complement();
void initialize_data_structures(data_structures *hashes);
memory_pool_collection* initialize_memory_pool_collection();
void initialize_statistics(statistics *stats);
void initialize_unit_sizes();
void initseq2Code();
void initdiff2hamming(unsigned char *difference);
int is_directory(const char *path);
void merge_feature_counts(gpointer key, gpointer value, gpointer user_data);
void merge_feature_sequences(gpointer key, gpointer value, gpointer user_data);
void merge_feature_umi_counts(gpointer key, gpointer value, gpointer user_data);
void merge_stats(statistics *merged_stats, statistics *thread_stats);
void merge_unmatched_barcodes(unmatched_barcodes_features_block_list *merged_list, unmatched_barcodes_features_block_list *thread_list);
void organize_fastq_files_by_directory(int positional_arg_count, int argc, char *argv[], int optind, char *barcodeFastqFilesString, char *forwardFastqFilesString, char *reverseFastqFilesString, fastq_files_collection *fastq_files, char *barcode_pattern, char *forward_pattern, char *reverse_pattern);
void organize_fastq_files_by_type(int positional_arg_count, int argc, char *argv[], int optind, char *barcodeFastqFilesString, char *forwardFastqFilesString, char *reverseFastqFilesString, fastq_files_collection *fastq_files, char *barcode_pattern, char *forward_pattern, char *reverse_pattern, int sample_flag);
void populate_sample_args(sample_args *args, int sample_index,char *directory, fastq_files_collection *fastq_files, feature_arrays *features, int maxHammingDistance, int nThreads, memory_pool_collection **pools, statistics *stats, data_structures *hashes, uint16_t stringency, uint16_t min_counts, int barcode_constant_offset, int feature_constant_offset, int read_buffer_lines, int average_read_length, double min_posterior, int consumer_threads_per_set);
void process_feature_sequence(char *sequence, feature_arrays *features, int maxHammingDistance, int nThreads, int feature_constant_offset, int max_feature_n, uint32_t *feature_index, int *hamming_distance, char *matching_sequence, uint16_t *match_position);
void process_files_in_sample(sample_args *args);
feature_arrays* read_features_file(const char* filename);
void reverse_complement_in_place(char *seq);
void reverse_complement_sequence(char *sequence,  char *reverse, int length);
void reverse_in_place(char *str);
void sort_samples_by_size(fastq_files_collection *fastq_files, int *sample_order);
int string2code(char *string, int sequence_length, unsigned char *code);
void update_feature_counts_from_code(unsigned char *code, char *umi, uint32_t feature_index, data_structures *hashes, memory_pool_collection *pools);
void update_umi_counts(unsigned char *code, char *umi,  uint32_t feature_index,data_structures *hashes, memory_pool_collection *pools);
void destroy_data_structures(data_structures *hashes);
void free_unmatched_barcodes_features_list(unmatched_barcodes_features_block_list *list);
unsigned char* read_whiteList(char *whitelist_filename,GHashTable *hash, int reverse_complement_flag);
int insert_feature_sequence(char *sequence, uint32_t feature_index, unsigned char hamming_distance, uint16_t match_position, data_structures *hashes, memory_pool_collection *pools);
int split_line(char *line, char *fields[], const char *split_string);

#endif // PROTOTYPES_H