// defines.h
#ifndef DEFINES_H
#define DEFINES_H

//defines that affect structs
#define BARCODE_LENGTH 16
#define UMI_LENGTH 12
#define BARCODE_CODE_LENGTH ((BARCODE_LENGTH + 3) / 4)
#define UMI_CODE_LENGTH ((UMI_LENGTH + 3) / 4)
#define MAX_FEATURES 128
#define MAX_BARCODE_MISMATCHES 3

//other defines
#define FEATURE_BARCODE_BLOCK_SIZE 30000
#define MAX_BARCODE_MISMATCHES 3
#define MIN_POSTERIOR 0.975
#define MAX_FEATURE_N 3  //customizable 
#define MAX_BARCODE_N 1  //customizable 
#define UMI_STORAGE_BLOCK 10000000
#define FEATURE_SEQUENCE_BLOCK_SIZE 100000
#define BARCODE_STORAGE_BLOCK_SIZE 30000 
#define MAX_MALLOCS 1000
#define MAX_FEATURE_SEQUENCE_LENGTH 64
#define READ_BUFFER_LINES 1024 //should be a multiple of 2 (we don't save first or 3rd line of fastq set of 4)
#define AVERAGE_READ_LENGTH 300 //make this generous since we will not be reallocing the circular read buffers


#define FILENAME_LENGTH 2048   
#define LINE_LENGTH 1024
#define NAME_LENGTH 40
#define MAX_FEATURE_CODE_LENGTH 40

#define FEATURE_NAME_LENGTH 32
#define LOOKUP_STRING "AAAAAAACAAAGAAATAACAAACCAACGAACTAAGAAAGCAAGGAAGTAATAAATCAATGAATTACAAACACACAGACATACCAACCCACCGACCTACGAACGCACGGACGTACTAACTCACTGACTTAGAAAGACAGAGAGATAGCAAGCCAGCGAGCTAGGAAGGCAGGGAGGTAGTAAGTCAGTGAGTTATAAATACATAGATATATCAATCCATCGATCTATGAATGCATGGATGTATTAATTCATTGATTTCAAACAACCAAGCAATCACACACCCACGCACTCAGACAGCCAGGCAGTCATACATCCATGCATTCCAACCACCCAGCCATCCCACCCCCCCGCCCTCCGACCGCCCGGCCGTCCTACCTCCCTGCCTTCGAACGACCGAGCGATCGCACGCCCGCGCGCTCGGACGGCCGGGCGGTCGTACGTCCGTGCGTTCTAACTACCTAGCTATCTCACTCCCTCGCTCTCTGACTGCCTGGCTGTCTTACTTCCTTGCTTTGAAAGAACGAAGGAATGACAGACCGACGGACTGAGAGAGCGAGGGAGTGATAGATCGATGGATTGCAAGCACGCAGGCATGCCAGCCCGCCGGCCTGCGAGCGCGCGGGCGTGCTAGCTCGCTGGCTTGGAAGGACGGAGGGATGGCAGGCCGGCGGGCTGGGAGGGCGGGGGGGTGGTAGGTCGGTGGGTTGTAAGTACGTAGGTATGTCAGTCCGTCGGTCTGTGAGTGCGTGGGTGTGTTAGTTCGTTGGTTTTAAATAACTAAGTAATTACATACCTACGTACTTAGATAGCTAGGTAGTTATATATCTATGTATTTCAATCACTCAGTCATTCCATCCCTCCGTCCTTCGATCGCTCGGTCGTTCTATCTCTCTGTCTTTGAATGACTGAGTGATTGCATGCCTGCGTGCTTGGATGGCTGGGTGGTTGTATGTCTGTGTGTTTTAATTACTTAGTTATTTCATTCCTTCGTTCTTTGATTGCTTGGTTGTTTTATTTCTTTGTTTT"
#define VISITED 255

// Debug print function

#define DEBUG_PRINT(fmt, ...) \
    do { \
        const char *env_debug = getenv("DEBUG"); \
        if ((env_debug && atoi(env_debug) != 0) || debug) { \
            fprintf(stderr, fmt, ##__VA_ARGS__); \
        } \
    } while (0)
#endif // DEFINES_H