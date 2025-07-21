#include "../include/io.h"

feature_arrays* read_features_file(const char* filename) {
    //expext a comma separated file with column names at least one with name and sequence fields
    int seq_size=0;
    int name_size=0;
    int code_size=0;
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open tags file");
        exit(EXIT_FAILURE);
    }
    char line[LINE_LENGTH];
    int count = 0;
    //skip the header
    //count the lines and check that the sequences are valid
    int maxFeatureLength=0;
    int seqIndex=-1;
    int nameIndex=-1;
    if (!fgets(line, LINE_LENGTH, file)) {
        perror("Failed to read tags header");
        exit(EXIT_FAILURE);
    }
    find_name_and_sequence_fields(line, &nameIndex, &seqIndex);
    
    while (fgets(line, LINE_LENGTH, file) != NULL) {
        get_feature_line_sizes(line, nameIndex, seqIndex, &name_size, &seq_size, &code_size, &maxFeatureLength);
        count++;
    }
    fprintf(stderr, "Read %d tags with max length %d\n", count, maxFeatureLength);
    feature_arrays *myfeatures = allocate_feature_arrays(name_size, seq_size, code_size, count, maxFeatureLength);
    //rewind the file and read the sequences
    fseek(file, 0, SEEK_SET);
    if (!fgets(line, LINE_LENGTH, file)) {
        perror("Failed to headers file");
        exit(EXIT_FAILURE);
    }
    count=0;
    while (fgets(line, LINE_LENGTH, file) != NULL) {
        process_feature_line(line, nameIndex, seqIndex, myfeatures, count);
        count++;
    }
    fprintf(stderr, "Read %d tags\n", count);
    fclose(file);
    return myfeatures;
}
void get_feature_line_sizes(char *line, int nameIndex, int seqIndex, int *name_size, int *seq_size, int *code_size, int *maxFeatureLength) {
    line[strcspn(line, "\r\n")] = '\0';
    char *fields[LINE_LENGTH];
    int nFields = split_line(line, fields, ",");
    if (seqIndex >= nFields || nameIndex >= nFields) {
        fprintf(stderr, "Error: Invalid line - too few fields %d \n",nFields);
        exit(EXIT_FAILURE);
    }
    char *tmpSeq = fields[seqIndex];
            // Remove possible newline character from the sequence
    if (!check_sequence(tmpSeq, strlen(tmpSeq))){
        //fprintf(stderr, "Error: Invalid sequence %s\n", tmpSeq);
        exit(EXIT_FAILURE);
    }
    const int string_length = strlen(tmpSeq);
    *seq_size += strlen(tmpSeq) + 1; 
    *code_size += string_length / 4;
    if (string_length % 4){
        (*code_size)++;
    } 
    if (strlen(tmpSeq) > *maxFeatureLength){
        *maxFeatureLength = strlen(tmpSeq);
    }
    *name_size += strlen(fields[nameIndex]) + 1;
}
void process_feature_line(char *line, int nameIndex, int seqIndex, feature_arrays *myfeatures, int count) {
    // Split the line by spaces and read the 3rd and 6th columns
    char *fields[LINE_LENGTH];
    line[strcspn(line, "\r\n")] = 0;
    int nFields = split_line(line, fields, ",");
    if (seqIndex >= nFields || nameIndex >= nFields) {
        fprintf(stderr, "Error: Invalid line - two few fields %d \n",nFields);
        exit(EXIT_FAILURE);
    }
    //copy the name and sequence into the feature arrays
    char *tmpName = fields[nameIndex];
    strcpy(myfeatures->feature_names[count], fields[nameIndex]);
    if (count + 1 < myfeatures->number_of_features) {
        myfeatures->feature_names[count + 1] = myfeatures->feature_names[count] + strlen(tmpName) + 1;
    }
    char *tmpSeq = fields[seqIndex];
    strcpy(myfeatures->feature_sequences[count], tmpSeq);
    myfeatures->feature_lengths[count] = strlen(tmpSeq);
    myfeatures->feature_code_lengths[count] = string2code(tmpSeq, strlen(tmpSeq), myfeatures->feature_codes[count]);
    GBytes *key = g_bytes_new_static(myfeatures->feature_codes[count], myfeatures->feature_code_lengths[count]);
    g_hash_table_insert(feature_code_hash, key, GUINT_TO_POINTER(count + 1));
    if (count + 1 < myfeatures->number_of_features) {
        myfeatures->feature_sequences[count + 1] = myfeatures->feature_sequences[count] + strlen(tmpSeq) + 1;
        myfeatures->feature_codes[count + 1] = myfeatures->feature_codes[count] + myfeatures->feature_code_lengths[count];
    }
}
feature_arrays* allocate_feature_arrays(int name_size, int seq_size, int code_size, int count, int maxFeatureLength) {
        feature_arrays *myfeatures = malloc(sizeof(feature_arrays));
        if (myfeatures == NULL) {
            fprintf(stderr, "Failed to allocate memory for feature arrays\n");
            exit(EXIT_FAILURE);
        }
        memset(myfeatures, 0, sizeof(feature_arrays));
        myfeatures->max_length = maxFeatureLength;
        myfeatures->feature_names_storage = malloc(name_size);
        myfeatures->feature_sequences_storage = malloc(seq_size);
        myfeatures->feature_codes_storage = malloc(code_size);
        myfeatures->feature_names = malloc(count * sizeof(char*));
        myfeatures->feature_lengths = malloc(count * sizeof(unsigned int));
        myfeatures->feature_code_lengths = malloc(count * sizeof(unsigned char));
        myfeatures->feature_codes = malloc(count * sizeof(unsigned char*));
        myfeatures->feature_sequences = malloc(count * sizeof(char*));
        myfeatures->number_of_features = count;

        // Check if any of the mallocs failed by checking for NULL pointers
        if (myfeatures->feature_names_storage == NULL || myfeatures->feature_sequences_storage == NULL || myfeatures->feature_codes_storage == NULL || myfeatures->feature_names == NULL || myfeatures->feature_lengths == NULL || myfeatures->feature_code_lengths == NULL || myfeatures->feature_codes == NULL) {
            fprintf(stderr, "Failed to allocate memory for feature arrays\n");
            exit(EXIT_FAILURE);
        }
        memset(myfeatures->feature_names_storage, 0, name_size);
        memset(myfeatures->feature_sequences_storage, 0, seq_size);
        memset(myfeatures->feature_codes_storage, 0, code_size);
        myfeatures->feature_names[0] = myfeatures->feature_names_storage;
        myfeatures->feature_sequences[0] = myfeatures->feature_sequences_storage;
        myfeatures->feature_codes[0] = myfeatures->feature_codes_storage;

        return myfeatures;
    }
void find_name_and_sequence_fields(char *line, int *nameIndex, int *seqIndex) {
    char *fields[LINE_LENGTH];
    //make sure that the line does not have a line feed
    line[strcspn(line, "\r\n")] = 0;
    int nFields = split_line(line, fields, ",");
    if (nFields < 2) {
        fprintf(stderr, "Error: Invalid header in tags file - there must be at least a name and sequence field \n");
        exit(EXIT_FAILURE);
    } else {
        for (int i = 0; i < nFields; i++) {
            if (strcmp(fields[i], "sequence") == 0) {
                *seqIndex = i;
            } else if (strcmp(fields[i], "name") == 0) {
                *nameIndex = i;
            }
            if (*seqIndex >= 0 && *nameIndex >= 0) {
                break;
            }
        }
    }
    if (*seqIndex < 0 || *nameIndex < 0) {
        fprintf(stderr, "Error: Invalid header in tags file - there must be at least a name and sequence field \n");
        exit(EXIT_FAILURE);
    }
}
int put_fastq_files_string_into_collection(char *fastqFilesString, char **fastq_files, int *nFiles, char *concatenated_fastq) {
    if (!fastqFilesString) {
        return 0;
    }
    char *this_file = concatenated_fastq;
    *nFiles = 1;
    strcpy(this_file, fastqFilesString);
    while (*this_file) {
        if (*this_file == ',') {
            (*nFiles)++;
            *this_file = '\0';
        }
        this_file++;
    }
    this_file = concatenated_fastq;
    for (int i = 0; i < *nFiles; i++) {
        fastq_files[i] = this_file;
        this_file += strlen(this_file) + 1;
    }
    return *nFiles;
}
void check_filecounts(fastq_files_collection *fastq_files) {
    if (!fastq_files->nbarcode_files) {
        fprintf(stderr, "Error: No barcode fastq files\n");
        exit(EXIT_FAILURE);
    }
    if (!fastq_files->nforward_files && !fastq_files->nreverse_files) {
        fprintf(stderr, "Error: No forward or reverse fastq files\n");
        exit(EXIT_FAILURE);
    }
    if (fastq_files->nforward_files && fastq_files->nforward_files != fastq_files->nbarcode_files) {
        fprintf(stderr, "Error: Unequal number of barcode and forward fastq files\n");
        exit(EXIT_FAILURE);
    }
    if (fastq_files->nreverse_files && fastq_files->nreverse_files != fastq_files->nbarcode_files) {
        fprintf(stderr, "Error: Unequal number of barcode and reverse fastq files\n");
        exit(EXIT_FAILURE);
    }
}
char* extract_sample_name(char *filename, char *pattern) {
    const char *samplename=get_basename(filename);
    char *pattern_position = strstr(samplename, pattern);
    if (!pattern_position) {
        fprintf(stderr, "Error: FASTQ file name %s does not contain the sample pattern %s\n", samplename, pattern);
        exit(EXIT_FAILURE);
    }
    *pattern_position = '\0';
    //find last underscore
    char *underscore_position = strrchr(samplename, '_');
    if (!underscore_position) {
        return (char*) samplename;
    }
    if (*(underscore_position+1) == 'L' && isdigit(*(underscore_position+2))){
        *underscore_position = '\0';
    }
    return (char*) samplename;
}
int count_character(char *string, char character) {
    int count = 0;
    for (int i = 0; string[i]; i++) {
        if (string[i] == character) {
            count++;
        }
    }
    return count;
}
int compare_filenames(const void *a, const void *b) {
    return strcmp(*(const char **)a, *(const char **)b);
}
int count_files_with_pattern(const char *directory_path, const char *pattern) {
    DIR *dir = opendir(directory_path);
    if (!dir) {
        perror("Unable to open directory");
        return -1;  // Return -1 to indicate an error
    }

    struct dirent *entry;
    struct stat file_stat;
    char filepath[FILENAME_LENGTH];
    int count = 0;

    // Iterate through the directory entries
    while ((entry = readdir(dir)) != NULL) {
        // Skip the "." and ".." entries
        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
            continue;
        }

        // Build the full file path
        snprintf(filepath, FILENAME_LENGTH, "%s/%s", directory_path, entry->d_name);

        // Check if the entry is a regular file (not a directory)
        if (stat(filepath, &file_stat) == 0 && S_ISREG(file_stat.st_mode)) {
            // Check if the filename contains the pattern
            if (strstr(entry->d_name, pattern) != NULL) {
                count++;  // Increment the count if a match is found
            }
        }
    }

    closedir(dir);
    return count;  // Return the count of matching files
}
char **find_files_with_pattern(const char *directory_path, const char *pattern, int *num_files_found) {
    DIR *dir = opendir(directory_path);
    if (!dir) {
        perror("Unable to open directory");
        return NULL;
    }

    struct dirent *entry;
    struct stat file_stat;
    char filepath[FILENAME_LENGTH];
    char **filepaths = NULL;  // Dynamically allocated array to store filepaths
    *num_files_found = 0;

    while ((entry = readdir(dir)) != NULL) {
        // Skip the "." and ".." entries
        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
            continue;
        }

        // Build the full file path
        strcpy(filepath, directory_path);
        if (filepath[strlen(filepath) - 1] != '/') {
            strcat(filepath, "/");
        }
        strcpy(filepath+strlen(filepath), entry->d_name);

        // Check if the entry is a regular file (not a directory)
        if (stat(filepath, &file_stat) == 0 && S_ISREG(file_stat.st_mode)) {
            // Check if the filename contains the pattern
            if (strstr(entry->d_name, pattern) != NULL) {
                // Store the matching filepath in the array
                filepaths = realloc(filepaths, (*num_files_found + 1) * sizeof(char *));
                if (!filepaths) {
                    perror("Memory allocation error");
                    exit(EXIT_FAILURE);
                }

                filepaths[*num_files_found] = strdup(filepath);  // Save a copy of the full file path
                if (!filepaths[*num_files_found]) {
                    perror("Memory allocation error");
                    exit(EXIT_FAILURE);
                }

                (*num_files_found)++;
            }
        }
    }

    closedir(dir);
    if (*num_files_found == 0) {
        free(filepaths);
        return NULL;  // Return NULL if no files were found
    }

    // Sort the filepaths alphabetically
    qsort(filepaths, *num_files_found, sizeof(char *), compare_filenames);

    return filepaths;  // Return the array of filepaths
}
void organize_fastq_files_by_directory(int positional_arg_count, int argc, char *argv[], int optind, char *barcodeFastqFilesString, char *forwardFastqFilesString, char *reverseFastqFilesString, fastq_files_collection *fastq_files, char *barcode_pattern, char *forward_pattern, char *reverse_pattern) {
    fastq_files->barcode_fastq = 0;
    fastq_files->forward_fastq = 0;
    fastq_files->reverse_fastq = 0;
    fastq_files->nbarcode_files = 0;
    fastq_files->nforward_files = 0;
    fastq_files->nreverse_files = 0;
    if (positional_arg_count) {
        //count the files in the first directory
        int barcode_file_exist=count_files_with_pattern(argv[optind], barcode_pattern);
        int forward_file_exist=count_files_with_pattern(argv[optind], forward_pattern);
        int reverse_file_exist=count_files_with_pattern(argv[optind], reverse_pattern);
        int total_barcode_files_found=0;
        for (int i=0;i < positional_arg_count;i++){
            total_barcode_files_found+=count_files_with_pattern(argv[optind+i], barcode_pattern);
        }
        if (!barcode_file_exist) {
            fprintf(stderr, "Error: No barcode fastq files found in directory %s\n", argv[optind]);
            exit(EXIT_FAILURE);
        }       
        if (!forward_file_exist && !reverse_file_exist) {
            fprintf(stderr, "Error: No forward or reverse fastq files found in directory %s\n", argv[optind]);
            exit(EXIT_FAILURE);
        }
        if (forward_file_exist) {
            fastq_files->forward_fastq = calloc(total_barcode_files_found, sizeof(char *));
            if (fastq_files->forward_fastq == NULL) {
                perror("Failed to allocate memory for forward fastq files");
                exit(EXIT_FAILURE);
            }       
        }
        if (reverse_file_exist) {
            fastq_files->reverse_fastq = calloc(total_barcode_files_found, sizeof(char *));
            if (fastq_files->reverse_fastq == NULL) {
                perror("Failed to allocate memory for reverse fastq files");
                exit(EXIT_FAILURE);
            }
        }
        fastq_files->barcode_fastq = calloc(total_barcode_files_found, sizeof(char *));
        fastq_files->sample_sizes=calloc(positional_arg_count,sizeof(int));
        fastq_files->sample_names=malloc(positional_arg_count*sizeof(char*));
        fastq_files->sample_offsets=malloc(positional_arg_count*sizeof(int));
        fastq_files->nsamples=positional_arg_count;
        fastq_files->sorted_index=malloc(positional_arg_count*sizeof(int));
        //check that the memory allocation was successful
        if (fastq_files->barcode_fastq == NULL || (forward_file_exist && fastq_files->forward_fastq == NULL) || (reverse_file_exist && fastq_files->reverse_fastq == NULL) || !fastq_files->sample_sizes || !fastq_files->sample_names || !fastq_files->sample_offsets || !fastq_files->sorted_index) {
            perror("Failed to allocate memory for fastq files");
            exit(EXIT_FAILURE);
        }
        total_barcode_files_found=0;
        for (int i = 0; i < positional_arg_count; i++) {
            int num_barcode_files_found = 0;
            int num_forward_files_found = 0;
            int num_reverse_files_found = 0;
            char *directory = strdup(argv[optind + i]); //this gets modified later
            char **sample_barcode_fastq=0, **sample_forward_fastq=0, **sample_reverse_fastq=0;

            sample_barcode_fastq=find_files_with_pattern(directory,barcode_pattern, &num_barcode_files_found);
            if (forward_file_exist) {
                sample_forward_fastq=find_files_with_pattern(directory,forward_pattern, &num_forward_files_found);
            }
            if (reverse_file_exist) {
                sample_reverse_fastq=find_files_with_pattern(directory,reverse_pattern, &num_reverse_files_found);
            }
            if (!num_barcode_files_found) {
                fprintf(stderr, "Error: No barcode fastq files found in directory %s\n", directory);
                exit(EXIT_FAILURE);
            }
            if (!num_forward_files_found && !num_reverse_files_found) {
                fprintf(stderr, "Error: No forward or reverse fastq files found in directory %s\n", directory);
                exit(EXIT_FAILURE);
            }
            if (num_forward_files_found && num_forward_files_found != num_barcode_files_found) {
                fprintf(stderr, "Error: Unequal number of barcode and forward fastq files in directory %s\n", directory);
                exit(EXIT_FAILURE);
            }
            if (num_reverse_files_found && num_reverse_files_found != num_barcode_files_found) {
                fprintf(stderr, "Error: Unequal number of barcode and reverse fastq files in directory %s\n", directory);
                exit(EXIT_FAILURE);
            }
            fastq_files->sample_sizes[i]=num_barcode_files_found;
            fastq_files->sample_offsets[i]=total_barcode_files_found;
            total_barcode_files_found+=num_barcode_files_found;
            for(int j=0;j<num_barcode_files_found;j++){
                fastq_files->barcode_fastq[fastq_files->nbarcode_files++]=sample_barcode_fastq[j];
                if (num_forward_files_found){
                    fastq_files->forward_fastq[fastq_files->nforward_files++]=sample_forward_fastq[j];
                }
                if (num_reverse_files_found){
                    fastq_files->reverse_fastq[fastq_files->nreverse_files++]=sample_reverse_fastq[j];
                }
            }
            
            fastq_files->sample_names[i] = strdup(get_basename(directory));
            fprintf(stderr, "directory %s Sample name %s\n", directory,fastq_files->sample_names[i]);
            free(directory);
        }
    }    
    sort_samples_by_size(fastq_files, fastq_files->sorted_index);
    int max_sample_size=0;
    for (int i=0; i<fastq_files->nsamples; i++){
        if (fastq_files->sample_sizes[i]>max_sample_size){
            max_sample_size=fastq_files->sample_sizes[i];
        }
    }
    fastq_files->max_sample_size=max_sample_size;
}

void organize_fastq_files_by_type(int positional_arg_count, int argc, char *argv[], int optind, char *barcodeFastqFilesString, char *forwardFastqFilesString, char *reverseFastqFilesString, fastq_files_collection *fastq_files, char *barcode_pattern, char *forward_pattern, char *reverse_pattern, int sample_flag) {
    fastq_files->barcode_fastq = 0;
    fastq_files->forward_fastq = 0;
    fastq_files->reverse_fastq = 0;
    fastq_files->nbarcode_files = 0;
    fastq_files->nforward_files = 0;
    fastq_files->nreverse_files = 0;
    if (positional_arg_count) {
        //allocate memory for the fastq files for now

        //count the size of the positional arguments
        size_t block_size = 0;
        for (int i = 0; i < positional_arg_count; i++) {
            block_size += strlen(argv[optind + i]) + 1;
        }
        fastq_files->concatenated_files = malloc(block_size);
        //check if the memory allocation was successful
        if (fastq_files->concatenated_files == NULL) {
            perror("Failed to allocate memory for concatenated fastq files");
            exit(EXIT_FAILURE);
        }

        //copy the positional arguments to the concatenated fastq files
        memset(fastq_files->concatenated_files, 0, block_size);
        char *this_file = fastq_files->concatenated_files;
        //count the file types and allocate as necessary
        for (int i = 0; i < positional_arg_count; i++) {
            strcpy(this_file, argv[optind + i]);
            char *barcode_pattern_position = strstr(this_file, barcode_pattern);
            char *forward_pattern_position = strstr(this_file, forward_pattern);
            char *reverse_pattern_position = strstr(this_file, reverse_pattern);
            if (!barcode_pattern_position && !forward_pattern_position && !reverse_pattern_position) {
                fprintf(stderr, "Error: FASTQ file name %s does not contain the barcode, forward and reverse patterns\n", this_file);
                exit(EXIT_FAILURE);
            }
            if (barcode_pattern_position && !forward_pattern_position && !reverse_pattern_position) {
                if (!fastq_files->barcode_fastq){
                    fastq_files->barcode_fastq = calloc(positional_arg_count, sizeof(char *));
                    if (fastq_files->barcode_fastq == NULL) {
                        perror("Failed to allocate memory for barcode fastq files");
                        exit(EXIT_FAILURE);
                    }
                }
                fastq_files->barcode_fastq[fastq_files->nbarcode_files++] = this_file;
            } else if (!barcode_pattern_position && forward_pattern_position && !reverse_pattern_position) {
                if (!fastq_files->forward_fastq){
                    fastq_files->forward_fastq = calloc(positional_arg_count, sizeof(char *));
                    if (fastq_files->forward_fastq == NULL) {
                        perror("Failed to allocate memory for forward fastq files");
                        exit(EXIT_FAILURE);
                    }
                }
                fastq_files->forward_fastq[fastq_files->nforward_files++] = this_file;
            } else if (!barcode_pattern_position && !forward_pattern_position && reverse_pattern_position) {
                if (!fastq_files->reverse_fastq){
                    fastq_files->reverse_fastq = calloc(positional_arg_count, sizeof(char *));
                    if (fastq_files->reverse_fastq == NULL) {
                        perror("Failed to allocate memory for reverse fastq files");
                        exit(EXIT_FAILURE);
                    }
                }
                fastq_files->reverse_fastq[fastq_files->nreverse_files++] = this_file;
            } else {
                fprintf(stderr, "Error: FASTQ file name %s contains more than one of the barcode, forward and reverse patterns\n", this_file);
                exit(EXIT_FAILURE);
            }
            this_file += strlen(argv[optind + i]) + 1;
        }
    } else {
        //find number of files in each type
        int barcode_count = 0;
        int forward_count = 0;
        int reverse_count = 0;
        if (barcodeFastqFilesString) {
            barcode_count = count_character(barcodeFastqFilesString, ',') + 1;
            fastq_files->barcode_fastq = calloc(barcode_count, sizeof(char *));
            if (fastq_files->barcode_fastq == NULL) {
                perror("Failed to allocate memory for barcode fastq files");
                exit(EXIT_FAILURE);
            }
        }
        if (forwardFastqFilesString) {
            forward_count = count_character(forwardFastqFilesString, ',') + 1;
            fastq_files->forward_fastq = calloc(forward_count, sizeof(char *));
            if (fastq_files->forward_fastq == NULL) {
                perror("Failed to allocate memory for forward fastq files");
                exit(EXIT_FAILURE);
            }
        }
        if (reverseFastqFilesString) {
            reverse_count = count_character(reverseFastqFilesString, ',') + 1;
            fastq_files->reverse_fastq = calloc(reverse_count, sizeof(char *));
            if (fastq_files->reverse_fastq == NULL) {
                perror("Failed to allocate memory for reverse fastq files");
                exit(EXIT_FAILURE);
            }
        }
        //find concatenated string length
        size_t fwd_len = forwardFastqFilesString ? strlen(forwardFastqFilesString) : 0;
        size_t rev_len = reverseFastqFilesString ? strlen(reverseFastqFilesString) : 0;
        size_t bar_len = barcodeFastqFilesString ? strlen(barcodeFastqFilesString) : 0;
        int concatenated_length = bar_len + fwd_len + rev_len + 3;
        fastq_files->concatenated_files = malloc(concatenated_length);
        if (fastq_files->concatenated_files == NULL) {
            perror("Failed to allocate memory for concatenated fastq files");
            exit(EXIT_FAILURE);
        }
        char *this_file = fastq_files->concatenated_files;
        if (barcodeFastqFilesString) {
            put_fastq_files_string_into_collection(barcodeFastqFilesString, fastq_files->barcode_fastq, &fastq_files->nbarcode_files, this_file);
        } else {
            perror("Must have some barcode fastq files");
            exit(EXIT_FAILURE);
        }
        put_fastq_files_string_into_collection(forwardFastqFilesString, fastq_files->forward_fastq, &fastq_files->nforward_files, this_file);
        put_fastq_files_string_into_collection(reverseFastqFilesString, fastq_files->reverse_fastq, &fastq_files->nreverse_files, this_file);
    }
    check_filecounts(fastq_files);
    fastq_files->nsamples=0;
    int name_length=0;
    for (int i = 0; i < fastq_files->nbarcode_files; i++) {
        name_length+=strlen(fastq_files->barcode_fastq[i])+1;
    }
    fprintf(stderr, "Name length %d\n", name_length);
    fastq_files->concatenated_sample_names=malloc(name_length+1);
    fastq_files->sample_sizes=malloc(fastq_files->nbarcode_files*sizeof(int));
    fastq_files->sample_names=malloc(fastq_files->nbarcode_files*sizeof(char*));
    fastq_files->sample_offsets=malloc(fastq_files->nbarcode_files*sizeof(int));
    //check if the memory allocation was successful
    if (fastq_files->concatenated_sample_names == NULL || fastq_files->sample_sizes == NULL) {
        perror("Failed to allocate memory for sample names");
        exit(EXIT_FAILURE);
    }
    memset(fastq_files->concatenated_sample_names, 0, strlen(fastq_files->concatenated_files)+1);
    char name_copy[name_length+1];
    char *last_sample=fastq_files->concatenated_sample_names;
    char *sample_name=0;
    fastq_files->nsamples=0;
    fastq_files->sample_offsets[0]=0;
    if (sample_flag){
        //assume that files are sorted by name
        for (int i = 0; i < fastq_files->nbarcode_files; i++) {
            strcpy(name_copy, fastq_files->barcode_fastq[i]); //copy is nccecessary because extract_sample_name modifies the string
            sample_name=extract_sample_name(name_copy, barcode_pattern);
            fprintf(stderr, "Sample name %s\n", sample_name);
            fprintf(stderr, "last sample name %s\n", last_sample);
            if (!fastq_files->nsamples){
                strcpy(last_sample, sample_name);
                fastq_files->sample_sizes[0]=1;
                fastq_files->sample_names[0]=last_sample;
                fastq_files->nsamples=1;
            }
            else if (strcmp(last_sample, sample_name)){
                last_sample+=strlen(last_sample)+1;
                fastq_files->sample_offsets[fastq_files->nsamples]=i;
                fastq_files->sample_names[fastq_files->nsamples]=last_sample;
                strcpy(last_sample, sample_name);
                fastq_files->sample_sizes[fastq_files->nsamples]=1;
                fastq_files->nsamples++;
            }
            else{
                fastq_files->sample_sizes[fastq_files->nsamples-1]++;
            }
        }
    }
    else{
        fastq_files->nsamples=1;
        fastq_files->sample_sizes[0]=fastq_files->nbarcode_files;
    }
    //now allocate the sample buffers
    fastq_files->sorted_index=malloc(fastq_files->nsamples*sizeof(int));
    if (!fastq_files->sorted_index){
        perror("Failed to allocate memory for sorted index");
        exit(EXIT_FAILURE);
    }
    sort_samples_by_size(fastq_files, fastq_files->sorted_index);

    int max_sample_size=0;
    for (int i=0; i<fastq_files->nsamples; i++){
        if (fastq_files->sample_sizes[i]>max_sample_size){
            max_sample_size=fastq_files->sample_sizes[i];
        }
    }
    fastq_files->max_sample_size=max_sample_size;
    //check that memory allocations are non null and free them
}

void populate_sample_args(sample_args *args, int sample_index,char *directory, fastq_files_collection *fastq_files, feature_arrays *features, int maxHammingDistance, int nThreads, memory_pool_collection **pools, statistics *stats, data_structures *hashes, uint16_t stringency, uint16_t min_counts, int barcode_constant_offset, int feature_constant_offset, int read_buffer_lines, int average_read_length, double min_posterior, int consumer_threads_per_set){ 
    args->sample_index = sample_index;
    args->directory = directory;
    args->fastq_files = fastq_files;
    args->features = features;
    args->maxHammingDistance = maxHammingDistance;
    args->nThreads = nThreads;
    args->pools = pools;
    args->stats = stats;
    args->hashes = hashes;
    args->stringency = stringency;
    args->min_counts = min_counts;
    args->barcode_constant_offset = barcode_constant_offset;
    args->feature_constant_offset = feature_constant_offset;
    args->average_read_length = average_read_length;
    args->read_buffer_lines = read_buffer_lines;
    args->min_posterior=min_posterior;
    args->consumer_threads_per_set=consumer_threads_per_set;

}

int find_number_of_fastq_files(int positional_arg_count,char *barcodeFastqFilesString, char *forwardFastqFilesString, char *reverseFastqFilesString){
    if (positional_arg_count){
        return positional_arg_count;
    }
    if (barcodeFastqFilesString == NULL){
        perror("No barcode fastq files specified");
        exit(EXIT_FAILURE);
    }
    int nFiles=1;
    if (forwardFastqFilesString != NULL){
        nFiles++;
    }
    if (reverseFastqFilesString != NULL){
        nFiles++;
    }
    for (int i=0; barcodeFastqFilesString[i]; i++){
        if (barcodeFastqFilesString[i] == ','){
            nFiles++;
        }
    }
    for (int i=0; forwardFastqFilesString[i]; i++){
        if (forwardFastqFilesString[i] == ','){
            nFiles++;
        }
    }
    for (int i=0; forwardFastqFilesString[i]; i++){
        if (forwardFastqFilesString[i] == ','){
            nFiles++;
        }
    }
    return nFiles;
}
int file_exists(const char *filename){
    struct stat buffer;
    return (stat(filename, &buffer) == 0);
}
const char* get_basename(const char* path) {
    DEBUG_PRINT( "Path %s\n", path);
    //case where the last character is a '/'
    if (path[0] && path[strlen(path)-1] == '/'){
        char *temp=strdup(path);
        temp[strlen(temp)-1]='\0'; //this should be cleaned at end but not a big deal
        return get_basename(temp);
    }
    const char *base = strrchr(path, '/');
    DEBUG_PRINT( "Base %s\n", base);
    return base ? base + 1 : path;  // If '/' is found, return the part after it. If not, return the original path.
}

void find_file_type(char* concatenated_patterns,int nPatterns){
    char patterns[nPatterns][strlen(concatenated_patterns)+1];
    char *this_pattern=concatenated_patterns;
    for (int i=0; i<nPatterns; i++){
        strcpy(patterns[i],this_pattern);
        this_pattern+=strlen(this_pattern)+1;
    }
    
}


