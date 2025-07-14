#ifndef IO_H
#define IO_H

#include "common.h"

// Function prototypes for file I/O operations
unsigned char* read_whiteList(char *whitelist_filename, GHashTable *hash, int reverse_complement_flag);
feature_arrays* read_features_file(const char* filename);
void organize_fastq_files_by_directory(int positional_arg_count, int argc, char *argv[], int optind, char *barcodeFastqFilesString, char *forwardFastqFilesString, char *reverseFastqFilesString, fastq_files_collection *fastq_files, char *barcode_pattern, char *forward_pattern, char *reverse_pattern);
void organize_fastq_files_by_type(int positional_arg_count, int argc, char *argv[], int optind, char *barcodeFastqFilesString, char *forwardFastqFilesString, char *reverseFastqFilesString, fastq_files_collection *fastq_files, char *barcode_pattern, char *forward_pattern, char *reverse_pattern, int sample_flag);
void read_fastq_chunks(gzFile barcode_fastqgz, gzFile forward_fastqgz, gzFile reverse_fastqgz, char *barcode_lines[], char *forward_lines[], char *reverse_lines[], char *buffer, char *done);
void open_fastq_files(const char *barcode_fastq, const char *forward_fastq, const char *reverse_fastq, gzFile *barcode_fastqgz, gzFile *forward_fastqgz, gzFile *reverse_fastqgz);
char **find_files_with_pattern(const char *directory_path, const char *pattern, int *num_files_found);
int count_files_with_pattern(const char *directory_path, const char *pattern);
size_t get_file_size(char *filepath);
void printFeatureCounts(feature_arrays *features, int *deduped_counts, int *barcoded_counts,int **coexpression_counts, int **coexpression_histograms, char *directory, data_structures *hashes, statistics *stats, uint16_t stringency, uint16_t min_counts);
int print_feature_sequences(feature_arrays *features, int *total_counts, char *directory, data_structures *hashes);
void generate_heatmap(const char *directory, feature_arrays *features, int **coexpression_histograms);
int mkdir_p(const char *path);
int file_exists(const char *filename);
const char* get_basename(const char* path);
char* extract_sample_name(char *filename, char *pattern);
void *read_fastqs_by_set(void *arg);
void *read_fastqs(void *arg);
void *consume_reads(void *arg);

#endif // IO_H