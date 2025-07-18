#ifndef IO_H
#define IO_H

#include "common.h"
#include "prototypes.h"
#include "utils.h"
#include "memory.h"
#include "globals.h"
// Function prototypes for file I/O operations
unsigned char* read_whiteList(char *whitelist_filename,GHashTable *hash, int reverse_complement_flag);
feature_arrays* read_features_file(const char* filename);
void get_feature_line_sizes(char *line, int nameIndex, int seqIndex, int *name_size, int *seq_size, int *code_size, int *maxFeatureLength);
void process_feature_line(char *line, int nameIndex, int seqIndex, feature_arrays *myfeatures, int count);
feature_arrays* allocate_feature_arrays(int name_size, int seq_size, int code_size, int count, int maxFeatureLength);
void find_name_and_sequence_fields(char *line, int *nameIndex, int *seqIndex);
int insert_feature_sequence(char *sequence, uint32_t feature_index, unsigned char hamming_distance, uint16_t match_position, data_structures *hashes, memory_pool_collection *pools);
int split_line(char *line, char *fields[], const char *split_string);
void free_feature_arrays(feature_arrays *features);
void free_memory_pool_collection(memory_pool_collection *pools);
void free_memory_pool_storage(memory_pool *pool);
storage_block* allocate_storage_block(size_t block_size);
memory_pool* initialize_storage_pool(size_t block_size, size_t blocks_per_pool);
void expand_memory_pool(memory_pool *pool);
int put_fastq_files_string_into_collection(char *fastqFilesString, char **fastq_files, int *nFiles, char *concatenated_fastq);
void check_filecounts(fastq_files_collection *fastq_files);
int count_character(char *string, char character);
int compare_filenames(const void *a, const void *b);
int count_files_with_pattern(const char *directory_path, const char *pattern);
void organize_fastq_files_by_directory(int positional_arg_count, int argc, char *argv[], int optind, char *barcodeFastqFilesString, char *forwardFastqFilesString, char *reverseFastqFilesString, fastq_files_collection *fastq_files, char *barcode_pattern, char *forward_pattern, char *reverse_pattern);
void organize_fastq_files_by_type(int positional_arg_count, int argc, char *argv[], int optind, char *barcodeFastqFilesString, char *forwardFastqFilesString, char *reverseFastqFilesString, fastq_files_collection *fastq_files, char *barcode_pattern, char *forward_pattern, char *reverse_pattern, int sample_flag);
void populate_sample_args(sample_args *args, int sample_index,char *directory, fastq_files_collection *fastq_files, feature_arrays *features, int maxHammingDistance, int nThreads, memory_pool_collection *pools, statistics *stats, data_structures *hashes, uint16_t stringency, uint16_t min_counts, int barcode_constant_offset, int feature_constant_offset, int read_buffer_lines, int average_read_length, int parallel_by_file, double min_posterior, int consumer_threads_per_set);
void sort_samples_by_size(fastq_files_collection *fastq_files, int *sample_order);
int find_number_of_fastq_files(int positional_arg_count,char *barcodeFastqFilesString, char *forwardFastqFilesString, char *reverseFastqFilesString);
int file_exists(const char *filename);
void find_file_type(char* concatenated_patterns,int nPatterns);
const char* get_basename(const char *path);

#endif // IO_H