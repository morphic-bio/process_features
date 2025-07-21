#include "../include/common.h"
#include "../include/globals.h"
#include "../include/prototypes.h"
#include "../include/utils.h"
#include "../include/io.h"

int main(int argc, char *argv[])
{   
    omp_set_nested(1);
    initseq2Code();
    initcode2seq();
    initdiff2hamming(diff2Hamming);
    initialize_complement();
    feature_code_hash = g_hash_table_new_full(g_bytes_hash, g_bytes_equal, (GDestroyNotify)g_bytes_unref, NULL);
    feature_arrays *features=0;
    int reverse_complement_whitelist=0;
    char *barcodeFastqFilesString=0;
    char *forwardFastqFilesString=0;
    char *reverseFastqFilesString=0;
    char barcode_pattern[LINE_LENGTH], forward_pattern[LINE_LENGTH], reverse_pattern[LINE_LENGTH];
    strcpy(barcode_pattern, "_R1_");
    strcpy(forward_pattern, "_R2_");
    strcpy(reverse_pattern, "_R3_");
    
    int maxHammingDistance=1;
    char directory[LINE_LENGTH];
    char whitelist_filename[4096];
    char sample_flag=1;
    char keep_existing=0;
    uint16_t stringency=1;
    uint16_t min_counts=1;
    int read_buffer_lines=READ_BUFFER_LINES;
    int average_read_length=AVERAGE_READ_LENGTH;
    int feature_constant_offset=0;
    int barcode_constant_offset=0;
    double min_posterior=MIN_POSTERIOR;

    int max_concurrent_processes=8;
    int consumer_threads_per_set=1;
    int search_threads_per_consumer=4;
    int set_consumer_threads_per_set=0;
    int set_search_threads_per_consumer=0;

    static struct option long_options[] = {
        {"whitelist", required_argument, 0, 'w'},
        {"featurelist", required_argument, 0, 'f'},
        {"maxHammingDistance", required_argument, 0, 'm'},
        {"feature_constant_offset", required_argument, 0, 'o'},
        {"barcode_constant_offset", required_argument, 0, 'B'},
        {"threads", required_argument, 0, 't'},
        {"search_threads", required_argument, 0, 'S'},
        {"stringency", required_argument, 0, 's'},
        {"min_counts", required_argument, 0, 'i'},
        {"umi_length",required_argument , 0, 'u'},
        {"barcode_length",required_argument , 0, 'b'},
        {"directory", required_argument, 0, 'd'},
        {"debug", no_argument, 0, 'v'},
        {"as_named", no_argument, 0, 'a'},
        {"reverse_complement_whitelist", no_argument, 0, 'r'},
        {"keep_existing", no_argument, 0, 'k'},
        {"read_buffer_lines", required_argument, 0, 'R'},
        {"average_read_length", required_argument, 0, 'L'},
        {"min_posterior", required_argument, 0, 'M'},
        {"consumer_threads_per_set", required_argument, 0, 'c'},
        {"barcode_fastqs", required_argument, 0, 0},
        {"forward_fastqs", required_argument, 0, 1},
        {"reverse_fastqs", required_argument, 0, 2},
        {"max_barcode_mismatches", required_argument, 0, 3},
        {"feature_n", required_argument, 0, 4},
        {"barcode_n", required_argument, 0, 5},
        {"barcode_fastq_pattern", required_argument, 0, 6},
        {"forward_fastq_pattern", required_argument, 0, 7},
        {"reverse_fastq_pattern", required_argument, 0, 8},
        {"max_reads", required_argument, 0, 9},
        {"limit_search", required_argument, 0, 10},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "w:b:f:m:s:S:i:t:T:o:d:u:c:vakrDB:R:L:M:", long_options, &option_index)) != -1) {
        switch (c) {
            case 'w': strcpy(whitelist_filename, optarg); break;
            case 'b': barcode_length=atoi(optarg); barcode_code_length=(barcode_length+3)/4; break;
            case 'f': features=read_features_file(optarg); number_of_features=features->number_of_features; maximum_feature_length=features->max_length; feature_code_length=(maximum_feature_length+3)/4; break;
            case 'm': maxHammingDistance=atoi(optarg); break;
            case 's': stringency=(uint16_t)atoi(optarg); break;
            case 'S': set_search_threads_per_consumer=atoi(optarg); break;
            case 'i': min_counts=(uint16_t)atoi(optarg); break;
            case 'd': strcpy(directory, optarg); if (directory[strlen(directory)-1] != '/'){ strcat(directory, "/"); } break;
            case 'o': feature_constant_offset=atoi(optarg); break;
            case 't': max_concurrent_processes=atoi(optarg); break;
            case 'u': umi_length=atoi(optarg); umi_code_length=(umi_length+3)/4; break;
            case 'c': set_consumer_threads_per_set=atoi(optarg); break;
            case 'v': debug = 1; break;
            case 'a': sample_flag=0; break;
            case 'k': keep_existing=1; break;
            case 'r': reverse_complement_whitelist=1; break;
            case 'B': barcode_constant_offset=atoi(optarg); break;
            case 'R': read_buffer_lines=atoi(optarg); break;
            case 'L': average_read_length=atoi(optarg); break;
            case 'M': min_posterior=(double)atof(optarg); break;
            case 0: barcodeFastqFilesString = strdup(optarg); break;
            case 1: forwardFastqFilesString = strdup(optarg); break;
            case 2: reverseFastqFilesString = strdup(optarg); break;
            case 3: max_barcode_mismatches=atoi(optarg); break;    
            case 4: max_feature_n=atoi(optarg); break;
            case 5: max_barcode_n=atoi(optarg); break;
            case 6: strcpy(barcode_pattern, optarg); break;
            case 7: strcpy(forward_pattern, optarg); break;
            case 8: strcpy(reverse_pattern, optarg); break;
            case 9: max_reads=atoll(optarg); break;
            case 10: limit_search = atoi(optarg); break;
            default: fprintf(stderr, "Usage: %s [options]\n", argv[0]); return 1;
        }
    }

    int positional_arg_count = argc - optind;
    if (positional_arg_count > 0 && (barcodeFastqFilesString != NULL || forwardFastqFilesString != NULL || reverseFastqFilesString != NULL)) {
        fprintf(stderr, "Error: Cannot specify both positional arguments and --barcode_fastqs\n");
        return 1;
    }
    if (strcmp(barcode_pattern, forward_pattern) == 0 || strcmp(barcode_pattern, reverse_pattern) == 0 || strcmp(forward_pattern, reverse_pattern) == 0){
        fprintf(stderr, "Error: Barcode, forward, and reverse patterns must be different\n");
        return 1;
    }
    fastq_files_collection fastq_files;
    memset(&fastq_files, 0, sizeof(fastq_files));
    if (positional_arg_count && is_directory(argv[optind])){
        organize_fastq_files_by_directory(positional_arg_count, argc, argv, optind, barcodeFastqFilesString, forwardFastqFilesString, reverseFastqFilesString, &fastq_files, barcode_pattern, forward_pattern, reverse_pattern);
    }
    else{
        organize_fastq_files_by_type(positional_arg_count, argc, argv,optind, barcodeFastqFilesString, forwardFastqFilesString, reverseFastqFilesString, &fastq_files, barcode_pattern, forward_pattern, reverse_pattern,sample_flag);
    }
    whitelist_hash = g_hash_table_new(hash_int32, equal_int32 );
    read_whiteList(whitelist_filename, whitelist_hash, reverse_complement_whitelist);
    initialize_unit_sizes();
    const int nSamples=fastq_files.nsamples;
    if (set_search_threads_per_consumer){
        search_threads_per_consumer=set_search_threads_per_consumer;
    }
    if (set_consumer_threads_per_set){
        consumer_threads_per_set=set_consumer_threads_per_set;
    } 
    fprintf(stderr, "Using %d consumer threads and %d search threads per consumer. Max concurrent processes %d\n", consumer_threads_per_set, search_threads_per_consumer, max_concurrent_processes);
    int concurrent_processes=0;
    atomic_int *thread_counter=mmap(NULL, sizeof(atomic_int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    if (thread_counter == MAP_FAILED){
        perror("Failed to allocate memory for thread counter");
        exit(EXIT_FAILURE);
    }
    atomic_init(thread_counter, 0);
    int threads_per_set = consumer_threads_per_set * (1 + search_threads_per_consumer);
    fprintf(stderr,"threads per set %d\n", threads_per_set);
    for (int index=0; index<nSamples; index++){
        const int i=fastq_files.sorted_index[index];
    
        if (concurrent_processes >= max_concurrent_processes){
            wait(NULL);
            concurrent_processes--;
        }
        pid_t pid=fork();
        if (pid< 0){
            perror("Failed to fork");
            exit(EXIT_FAILURE); 
        }
        else if (pid== 0){
            atomic_fetch_add(thread_counter, threads_per_set);
            
            char sample_directory[FILENAME_LENGTH];
            if (sample_flag){
                strcpy(sample_directory, directory);
                strcat(sample_directory, fastq_files.sample_names[i]);
                strcat(sample_directory, "/");
            }
            else{
                strcpy(sample_directory, directory);
                strcat(sample_directory, fastq_files.sample_names[0]);
                strcat(sample_directory, "/");
            }
            fprintf(stderr, "Processing sample directory %s\n", sample_directory);
            if (existing_output_skip(keep_existing, sample_directory)) exit(0);
            sample_args args;
            // The following are now initialized in process_files_in_sample
            // memory_pool_collection *pools=initialize_memory_pool_collection();
            // statistics stats;
            // data_structures hashes;
            // initialize_data_structures(&hashes);
            // initialize_statistics(&stats);
            populate_sample_args(&args,i, sample_directory, &fastq_files, features, maxHammingDistance, search_threads_per_consumer, NULL, NULL, NULL, stringency, min_counts, barcode_constant_offset, feature_constant_offset, read_buffer_lines, average_read_length,min_posterior,consumer_threads_per_set);
            process_files_in_sample(&args);
            // cleanup_sample is handled within process_files_in_sample
            atomic_fetch_add(thread_counter, -threads_per_set);
            exit(0);
        }
        concurrent_processes++;

    }
    while (concurrent_processes > 0) {
        wait(NULL);
        concurrent_processes--;
    }
    g_hash_table_destroy(whitelist_hash);
    free(whitelist);
    if (barcodeFastqFilesString) free(barcodeFastqFilesString);
    if (forwardFastqFilesString) free(forwardFastqFilesString);
    if (reverseFastqFilesString) free(reverseFastqFilesString);
    free_feature_arrays(features);
    free_fastq_files_collection(&fastq_files);
    return 0;
}