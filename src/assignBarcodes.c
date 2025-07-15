#include "../include/common.h"
#include "../include/globals.h"
#include "../include/prototypes.h"
#include "../include/utils.h"
#include "../include/memory.h"

//will  print if DEBUG is set or debug=1
//code for feature sequences stats

//code for heatmap generation
#ifndef NO_HEATMAP
#include <cairo.h>
#include "../include/plasma_colormap_16.h"
#include "../include/plasma_colormap_64.h"
#include "../include/plasma_colormap_256.h"
#include "../include/plasma_colormap_1024.h"

const double (*select_colormap(int dynamic_range, int *colormap_size))[3] {
    if (dynamic_range < 20) {
        *colormap_size = 16;
        return plasma_colormap_16;
    } else if (dynamic_range < 100) {
        *colormap_size = 64;
        return plasma_colormap_64;
    } else if (dynamic_range < 1000) {
        *colormap_size = 256;
        return plasma_colormap_256;
    } else {
        *colormap_size = 1024;
        return plasma_colormap_1024;
    }
}

// Function to map a normalized value (0–1) to an RGB color using a colormap
void value_to_color(double normalized_value, double *r, double *g, double *b, const double colormap[][3], int colormap_size) {
    if (normalized_value < 0.0) normalized_value = 0.0;
    if (normalized_value > 1.0) normalized_value = 1.0;
    int index = (int)(normalized_value * (colormap_size - 1));
    *r = colormap[index][0];
    *g = colormap[index][1];
    *b = colormap[index][2];
}

// Function to normalize values for coloring
double normalize(int value, int max_value) {
    return (double)value / max_value;
}

// Function to generate heatmap with dynamic padding
void generate_heatmap(const char *directory, feature_arrays *features, int **coexpression_histograms) {
    char output_file[1024];
    if (directory[strlen(directory) - 1] == '/') {
        snprintf(output_file, sizeof(output_file), "%sheatmap.png", directory);
    } else {
        snprintf(output_file, sizeof(output_file), "%s/heatmap.png", directory);
    }
    int max_name_length=0;
    const int num_rows = features->number_of_features;
    int column_sums[num_rows];
    memset(column_sums,0,num_rows*sizeof(int));
    for (int i=0; i<num_rows; i++){
        if (strlen(features->feature_names[i])>max_name_length){
            max_name_length=strlen(features->feature_names[i]);
        }
    }

    int left_padding = BASE_PADDING + max_name_length * 7;
    int right_padding = BASE_PADDING + 50;
    int filtered_rows = 0;
    int filter_mask[num_rows];
    int num_cols = 0;
    
    for (int i = 0; i < num_rows; i++) {
       filter_mask[i]=coexpression_histograms[i+1][0]>0;
         if (filter_mask[i]){
              filtered_rows++;
         }
       if (coexpression_histograms[i+1][0]> num_cols){
           num_cols=coexpression_histograms[i+1][0];
       }
    }
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            column_sums[j] += coexpression_histograms[i+1][j+1];
        }
    }
    for (int j = 0; j < num_cols; j++) {
        column_sums[j] /= j+1;
    }
    int max_column_sum = 0;
    for (int j = 0; j < num_cols; j++) {
        if (column_sums[j] > max_column_sum) {
            max_column_sum = column_sums[j];
        }
    }
    int width = left_padding + num_cols * CELL_SIZE + BAR_WIDTH + right_padding;
    int height = BASE_PADDING + BAR_GRAPH_HEIGHT + BASE_PADDING + 20 + filtered_rows * CELL_SIZE + BASE_PADDING;

    cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
    cairo_t *cr = cairo_create(surface);

    for (int j = 0; j < num_cols; j++) {
        double bar_height = (double)column_sums[j] / max_column_sum * BAR_GRAPH_HEIGHT;
        cairo_set_source_rgb(cr, 0.2, 0.4, 0.8);
        cairo_rectangle(cr, left_padding + j * CELL_SIZE, BASE_PADDING + BAR_GRAPH_HEIGHT - bar_height, CELL_SIZE, bar_height);
        cairo_fill(cr);
    }
    cairo_set_font_size(cr, 10);
    cairo_set_source_rgb(cr, 0, 0, 0);

    char max_label[20];
    snprintf(max_label, sizeof(max_label), "%-d", max_column_sum);
    cairo_move_to(cr, left_padding - 30 - 10, BASE_PADDING + 5);
    cairo_show_text(cr, max_label);

    char midpoint_label[20];
    int midpoint_value = max_column_sum / 2;
    snprintf(midpoint_label, sizeof(midpoint_label), "%d", midpoint_value);
    cairo_move_to(cr, left_padding - 30 - 10, BASE_PADDING + BAR_GRAPH_HEIGHT / 2 + 5);
    cairo_show_text(cr, midpoint_label);
    
    cairo_set_font_size(cr, 10);
    cairo_set_source_rgb(cr, 0, 0, 0);
    for (int j = 0; j < num_cols; j++) {
        if ((j+1) % 5 == 0) {
            char label[10];
            snprintf(label, sizeof(label), "%d", j+1);
            cairo_move_to(cr, left_padding + j * CELL_SIZE + CELL_SIZE / 4, BASE_PADDING + BAR_GRAPH_HEIGHT + 15);
            cairo_show_text(cr, label);
        }
    }

    int max_value = 0;
    for (int i = 0; i < num_rows; i++) {
        if (!filter_mask[i]) continue;
        for (int j = 0; j < num_cols; j++) {
            if (coexpression_histograms[i+1][j+1] > max_value) {
                max_value = coexpression_histograms[i+1][j+1];
            }
        }
    }
    int draw_row = 0;
    const double (*plasma_colormap)[3];
    int colormap_size;
    plasma_colormap = select_colormap(max_value, &colormap_size);

    for (int i = 0; i < num_rows; i++) {
        if (!filter_mask[i]) continue;
        for (int j = 0; j < num_cols; j++) {
            double intensity = normalize(coexpression_histograms[i+1][j+1], max_value);
            double r, g, b;
            value_to_color(intensity, &r, &g, &b, plasma_colormap, colormap_size);
            cairo_set_source_rgb(cr, r, g, b);
            cairo_rectangle(cr, left_padding + j * CELL_SIZE, BASE_PADDING + BAR_GRAPH_HEIGHT + 20 + draw_row * CELL_SIZE, CELL_SIZE, CELL_SIZE);
            cairo_fill(cr);
        }
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 10);
        cairo_move_to(cr, BASE_PADDING, BASE_PADDING + BAR_GRAPH_HEIGHT + 20 + draw_row * CELL_SIZE + CELL_SIZE / 2);
        cairo_show_text(cr, features->feature_names[i]);
        draw_row++;
    }
    for (int i = 0; i < height - BAR_GRAPH_HEIGHT - 20 - BASE_PADDING*3 ; i++) {
        double normalized_value = 1.0 - (double)i / (height - BAR_GRAPH_HEIGHT - BASE_PADDING * 3 - 20);
        double r, g, b;
        value_to_color(normalized_value, &r, &g, &b, plasma_colormap, colormap_size);
        cairo_set_source_rgb(cr, r, g, b);
        cairo_rectangle(cr, width - right_padding - BAR_WIDTH, BASE_PADDING + BAR_GRAPH_HEIGHT + 20 + i, BAR_WIDTH, 1);
        cairo_fill(cr);
    }

    cairo_set_font_size(cr, 10);
    for (int i = 0; i <= height - BAR_GRAPH_HEIGHT - BASE_PADDING * 3 - 20; i += 10 * CELL_SIZE) {
        double normalized_value = 1.0 - (double)i / (height - BAR_GRAPH_HEIGHT - BASE_PADDING * 3 - 20 );
        int value = (int)(normalized_value * max_value);
        char label[20];
        snprintf(label, sizeof(label), "%d", value);
        cairo_move_to(cr, width - right_padding + 5, BASE_PADDING + BAR_GRAPH_HEIGHT + 20 + i + 4);
        cairo_show_text(cr, label);
    }

    cairo_surface_write_to_png(surface, output_file);

    cairo_destroy(cr);
    cairo_surface_destroy(surface);
}
#endif

void destroy_feature_umi_counts(gpointer data) {
    feature_umi_counts *umi_counts = (feature_umi_counts*)data;
    if (umi_counts && umi_counts->counts) {
        g_hash_table_destroy(umi_counts->counts);
    }
    // We do NOT free(umi_counts) itself, because it was allocated
    // from our custom memory pool, which is freed all at once later.
}

void initialize_complement(){
    match['A']='T';
    match['T']='A';
    match['C']='G';
    match['G']='C';
    match['N']='N';
}
void initialize_statistics(statistics *stats) {
    stats->start_time = get_time_in_seconds();
    stats->nMismatches = 0;
    stats->recovered = 0;
    stats->pending = 0;
    stats->valid = 0;
    stats->pending_recovered = 0;
    stats->total_unmatched_features = 0;
    stats->number_of_reads = 0;
    stats->unmatched_list.first_entry = NULL;
    stats->unmatched_list.last_entry = NULL;
}

void free_unmatched_barcodes_features_list(unmatched_barcodes_features_block_list *list) {
    unmatched_barcodes_features_block *current = list->first_entry;
    while (current != NULL) {
        unmatched_barcodes_features_block *next = current->next;
        free(current);
        current = next;
    }
    list->first_entry = NULL;
    list->last_entry = NULL;
}
double get_time_in_seconds() {
    struct timeval time;
    gettimeofday(&time, NULL);
    return time.tv_sec + (time.tv_usec / 1000000.0);
}
int mkdir_p(const char *path) {
    char temp[1024];
    char *p = NULL;
    size_t len;

    // Copy path and ensure it ends with '/'
    snprintf(temp, sizeof(temp), "%s", path);
    len = strlen(temp);
    if (temp[len - 1] == '/')
        temp[len - 1] = 0;

    // Iterate through each directory in the path
    for (p = temp + 1; *p; p++) {
        if (*p == '/') {
            *p = 0;

            // Create directory if it doesn't exist
            if (mkdir(temp, S_IRWXU) != 0 && errno != EEXIST) {
                perror("mkdir");
                return -1;
            }
            *p = '/';
        }
    }
    // Create the final directory
    if (mkdir(temp, S_IRWXU) != 0 && errno != EEXIST) {
        perror("mkdir");
        return -1;
    }
    return 0;
}

void lowerCaseDifferences(char *ref, char *query, int length){
    for (int i=0; i<length; i++){
        if (ref[i] != query[i]){
            if (query[i] == '\n'){
               //put x until the length of the query and then add a \0
               for (int j=i; j<length; j++){
                   query[j]='x';
               }
               query[length]='\0';
               return;
            }
            query[i]=tolower(query[i]);
        }
    } 
}   
int compare_feature_sequences(const void *a, const void *b) {
    const feature_sequences *fa = *(const feature_sequences **)a;
    const feature_sequences *fb = *(const feature_sequences **)b;

    // First sort by feature_index
    if (fa->feature_index != fb->feature_index)
        return (fa->feature_index > fb->feature_index) - (fa->feature_index < fb->feature_index);

    // Then sort by hamming_distance
    if (fa->hamming_distance != fb->hamming_distance)
        return fa->hamming_distance - fb->hamming_distance;

    // Lastly sort by counts
    return fb->counts - fa->counts;
}

void initseq2Code(){
    for (int i=0; i<256; i++){
        seq2code[i]=0;
    }
    seq2code['A']=0;
    seq2code['C']=1;
    seq2code['G']=2;
    seq2code['T']=3;
}
void initcode2seq(){
    char temp[1025];
    strcpy(temp, LOOKUP_STRING);
    for (int i=0; i<256; i++){
        code2seq[i][0]=temp[4*i];
        code2seq[i][1]=temp[4*i+1];
        code2seq[i][2]=temp[4*i+2];
        code2seq[i][3]=temp[4*i+3];
    }
}
void initdiff2hamming(unsigned char *difference){
    memset(difference,0,256);
    for (int i=0; i<256; i++){
        difference[i]=0;
          unsigned char mask=3;
        //check if first 2 bits are 0
        for (int j=0; j<4; j++){
            if (i & mask){
                difference[i]++;
            }
            mask=mask << 2;
        }
        //DEBUG_PRINT( "Difference %d %d\n", i, difference[i]);
    }
}


void free_feature_arrays(feature_arrays *features) {
    free(features->feature_names_storage);
    free(features->feature_lengths);
    free(features->feature_code_lengths);
    free(features->feature_sequences_storage);
    free(features->feature_codes_storage);
    free(features->feature_names);
    free(features->feature_sequences);
    free(features->feature_codes);
    free(features);
}

void initialize_unit_sizes(){
    // Size of feature_counts (rounded up to alignment of uint16_t)
    size_t feature_counts_size = sizeof(feature_counts);
    size_t feature_counts_alignment = __alignof__(feature_counts);
    dynamic_struct_sizes.feature_counts = (feature_counts_size + feature_counts_alignment - 1) & ~(feature_counts_alignment - 1);


    // Size of feature_umi_counts (rounded up to alignment of uint16_t)
    size_t feature_umi_counts_size = sizeof(feature_umi_counts);
    size_t feature_umi_counts_alignment = __alignof__(feature_umi_counts);
    dynamic_struct_sizes.feature_umi_counts = (feature_umi_counts_size + feature_umi_counts_alignment - 1) & ~(feature_umi_counts_alignment - 1);


    // Size of feature_sequences (rounded up to alignment of char)
    size_t feature_sequences_size = sizeof(feature_sequences) + maximum_feature_length + 1;
    DEBUG_PRINT( "Feature sequences size %ld\n", feature_sequences_size);
    size_t feature_sequences_alignment = __alignof__(char);  // Use __alignof__ instead of alignof
    dynamic_struct_sizes.feature_sequences = (feature_sequences_size + feature_sequences_alignment - 1) & ~(feature_sequences_alignment - 1);
    DEBUG_PRINT( "Adjusted feature sequences size %ld\n", dynamic_struct_sizes.feature_sequences);

    // Size of unmatched_barcodes_features_block (rounded up to alignment of unsigned char)
    size_t unmatched_barcodes_features_block_size = sizeof(unmatched_barcodes_features_block)
            + sizeof(uint32_t) // For feature_index
            + sizeof(uint16_t) // For match_position
            + barcode_code_length 
            + umi_code_length
            + 1 // For number_of_closest_barcodes
            + (max_barcode_mismatches + 1) * barcode_code_length
            + (max_barcode_mismatches + 1);
    size_t unmatched_barcodes_features_block_alignment = __alignof__(uint32_t);  // Use __alignof__ instead of alignof
    dynamic_struct_sizes.unmatched_barcodes_features_block = (unmatched_barcodes_features_block_size + unmatched_barcodes_features_block_alignment - 1) & ~(unmatched_barcodes_features_block_alignment - 1);
}
int is_directory(const char *path) {
    struct stat path_stat;
    
    // Check if stat call is successful
    if (stat(path, &path_stat) != 0) {
        // Error: path not found or other issue
        perror("Error accessing path");
        return -1;
    }
    
    // Check if it is a directory
    if (S_ISDIR(path_stat.st_mode)) {
        return 1;
    }

    // Check if it is a file
    if (S_ISREG(path_stat.st_mode)) {
        return 0;
    }

    // In case it's neither a file nor a directory
    return -1;
}
void read_unmatched_features_block(unmatched_barcodes_features_block *entry_block, unmatched_barcodes_features *entry){
    entry->next=entry_block->next;
    unsigned char *storage = entry_block->storage;
    memcpy(&entry->feature_index, storage, sizeof(uint32_t));
    storage += sizeof(uint32_t);
    memcpy(&entry->match_position, storage, sizeof(uint16_t));
    storage += sizeof(uint16_t);
    entry->barcode=storage;
    entry->umi=entry->barcode+barcode_code_length;
    entry->number_of_closest_barcodes=entry->umi[umi_code_length];
    entry->closest_barcodes=entry->umi+umi_code_length+1;
    entry->Qscores=entry->closest_barcodes+(barcode_code_length)*(max_barcode_mismatches+1);
}
int insert_feature_sequence(char *sequence, uint32_t feature_index, unsigned char hamming_distance, uint16_t match_position, data_structures *hashes, memory_pool_collection *pools){
    gpointer *value=g_hash_table_lookup(hashes->unique_features_match, sequence);
    if (value){
        feature_sequences *entry = (feature_sequences*)value;
        entry->counts++;
        return 0;
    }
    else{
        feature_sequences *new_entry = (feature_sequences*) allocate_memory_from_pool(pools->feature_sequences_pool);
        strcpy(new_entry->sequence, sequence);
        new_entry->feature_index=feature_index;
        new_entry->hamming_distance=hamming_distance;
        new_entry->match_position=match_position;
        new_entry->counts=1;
        g_hash_table_insert(hashes->unique_features_match, new_entry->sequence, new_entry);
        return 1;
    }   
}
int print_feature_sequences(feature_arrays *features, int *total_counts, char *directory, data_structures *hashes){
    //remember that 1 is the first feature
    memset(total_counts, 0, features->number_of_features * sizeof(int));
    char feature_sequences_filename[FILENAME_LENGTH];
    sprintf(feature_sequences_filename, "%s/feature_sequences.txt", directory);
    FILE *feature_sequencesfp = fopen(feature_sequences_filename, "w");
    if (feature_sequencesfp == NULL) {
        fprintf(stderr, "Failed to open feature sequences file\n");
        exit(EXIT_FAILURE);
    }   
    GHashTableIter iter;
    gpointer key, value;
    int i=0;
    g_hash_table_iter_init(&iter, hashes->unique_features_match);
    int count = g_hash_table_size(hashes->unique_features_match);
    gpointer *array = g_new(gpointer, count);
    if(!array){
        fprintf(stderr, "Failed to allocate memory for array\n");
        exit(EXIT_FAILURE);
    }   
    while (g_hash_table_iter_next(&iter, &key, &value)) { 
        array[i++] = value;
    }
    qsort(array, count, sizeof(gpointer), compare_feature_sequences);
    fprintf(feature_sequencesfp, "Feature Index Sequence Hamming Distance Counts Match Position Feature Name\n");   
    for (i = 0; i < count; i++) {
        char annotated_sequence[LINE_LENGTH];
        const int mapped_index=((feature_sequences*)array[i])->feature_index-1;
        feature_sequences *entry = (feature_sequences*)array[i];
        total_counts[mapped_index]+=entry->counts;
        strcpy(annotated_sequence, entry->sequence);
        int feature_length=features->feature_lengths[entry->feature_index-1];
        lowerCaseDifferences(features->feature_sequences[entry->feature_index-1],annotated_sequence,feature_length);
        fprintf(feature_sequencesfp, "%7u %s %2d %7d %5u %s\n", entry->feature_index,annotated_sequence, entry->hamming_distance,entry->counts, entry->match_position, features->feature_names[entry->feature_index-1]);
    }
    g_free(array);
    return 0;
}



unmatched_barcodes_features_block* add_unmatched_barcode_store_feature(unsigned char *barcodes, unsigned char* corrected_barcodes, char *umi, unsigned char *qscores, uint32_t feature_index, int number_of_variants, uint16_t match_position, memory_pool_collection *pools, statistics *stats){
    unmatched_barcodes_features_block *new_entry = (unmatched_barcodes_features_block*)allocate_memory_from_pool(pools->unmatched_barcodes_features_block_pool); 

    if (stats->unmatched_list.first_entry== NULL){
        stats->unmatched_list.first_entry=new_entry;
        stats->unmatched_list.last_entry=new_entry;
    }
    else{
        stats->unmatched_list.last_entry->next=new_entry;
        stats->unmatched_list.last_entry=new_entry;
    }
    unsigned char *storage = new_entry->storage; 
    memcpy(storage, &feature_index, sizeof(uint32_t));
    storage += sizeof(uint32_t);
    memcpy(storage, &match_position, sizeof(uint16_t));
    storage += sizeof(uint16_t);
    memcpy(storage, barcodes, barcode_code_length);
    storage+=barcode_code_length;
    string2code(umi, umi_length, storage);
    storage+=umi_code_length;    
    storage[0]=number_of_variants;
    storage++;
    //leave empty space for the closest barcodes    
    unsigned char* qscores_storage = storage + (barcode_code_length)*(max_barcode_mismatches+1);
    memcpy(storage, corrected_barcodes, number_of_variants * barcode_code_length);
    memcpy(qscores_storage, qscores, number_of_variants);
    new_entry->next=NULL;
    return new_entry;
}


unsigned char* read_whiteList(char *whitelist_filename,GHashTable *hash, int reverse_complement_flag){
    FILE *whitelist_file = fopen(whitelist_filename, "r");
    if (whitelist_file == NULL) {
        fprintf(stderr, "Failed to open whitelist file %s\n", whitelist_filename);
        exit(EXIT_FAILURE);
    }
    char line[LINE_LENGTH];
    //just in case the file is corrupted - calculate the sizes by counting lines
    //find the number of lines in the file by counting the number of newlines
    char ch=0,lastChar=0;
    size_t line_count=0;
    //from the first line calculate the barcode length and check if the sequence is valid
    int sequence_count=0;
    while((ch = fgetc(whitelist_file)) != EOF) {
        if(ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'N'){
            sequence_count++;
        }
        else if (ch == '\n') {
            if(sequence_count){
                if (barcode_length==0){
                    barcode_length=sequence_count;
                    barcode_code_length=(barcode_length+3)/4;
                }
                else if (sequence_count != barcode_length){
                    fprintf(stderr, "Error: Invalid barcode length %d\n", sequence_count);
                    exit(EXIT_FAILURE);
                }
                sequence_count=0;
                line_count++;
            }
        }
        else{
            fprintf(stderr, "Error: Invalid character %c in whitelist file\n", ch);
            exit(EXIT_FAILURE);
        }
    }
    // If the last character isn't a newline, increment the count
    if (lastChar != '\n') {
        line_count++;
    }
    
    whitelist = malloc(line_count * barcode_code_length);
    if (whitelist == NULL) {
        fprintf(stderr, "Failed to allocate memory for whitelist storage\n");
        exit(EXIT_FAILURE);
    }
    memset(whitelist, 0, line_count * barcode_code_length);
    //reset the file pointer to the beginning of the file
    fseek(whitelist_file, 0, SEEK_SET);
    size_t nBarcodes=0;
    //set line to zero
    memset(line, 0, LINE_LENGTH);
    char barcode_sequence[LINE_LENGTH];
    while ( fgets(line, LINE_LENGTH, whitelist_file) != NULL) {
        if (!check_sequence(line, barcode_length)){
            fprintf(stderr, "Error: Invalid barcode sequence of expected length %d %s\n",barcode_length,line);
            exit(EXIT_FAILURE);
        }
        //up to 16 characters will be stored in uint32_t
        memset(line+barcode_length, 0, LINE_LENGTH-barcode_length);
        if (reverse_complement_flag){
            reverse_complement_sequence(line, barcode_sequence, barcode_length);
        }
        else{
            memcpy(barcode_sequence, line, barcode_length+1);
        }
        int j=0;
        int i=0;
        const size_t offset=nBarcodes*barcode_code_length;
        memset(whitelist+offset, 0, barcode_code_length);
        while(j<barcode_length){
            unsigned char char_value=seq2code[(unsigned char)barcode_sequence[j]]<<6 | seq2code[(unsigned char)barcode_sequence[j+1]]<<4 | seq2code[(unsigned char)barcode_sequence[j+2]]<<2 | seq2code[(unsigned char)barcode_sequence[j+3]];
            whitelist[offset+i]=char_value;
            i++;
            j+=4;
        }
        if(!g_hash_table_insert(hash,  (uint32_t*)(whitelist+offset) , whitelist+offset)){
            fprintf(stderr, "Error: Failed to insert barcode %s into the hash table %ld\n", line, nBarcodes);
            exit(EXIT_FAILURE); 
        }
        nBarcodes++;
    }
    fprintf(stderr, "Read %ld barcodes\n", nBarcodes);
    fclose(whitelist_file);
    return whitelist;
}
int split_line(char *line, char *fields[], const char *split_string) {
    int count = 0;
    char *token;

    // Use strtok to split the line by split_string
    token = strtok(line, split_string);
    while (token != NULL) {
        fields[count++] = token;  // Store pointer to the split string
        token = strtok(NULL, split_string);  // Get next token
    }
    return count;  // Return the number of fields
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
char check_sequence(char *sequence, int sequence_length){
    for (int i=0; i<sequence_length; i++){
        if (sequence[i] == 'A' || sequence[i] == 'C' || sequence[i] == 'G' || sequence[i] == 'T'){
            continue;
        }
        return 0;
    }
    return 1;
}
int string2code(char *string, int sequence_length, unsigned char *code){
    const int last_element=sequence_length/4;
    for (int i=0; i<last_element; i++){
        code[i]=seq2code[(unsigned char)string[4*i]]<<6 | seq2code[(unsigned char)string[4*i+1]]<<4 | seq2code[(unsigned char)string[4*i+2]]<<2 | seq2code[(unsigned char)string[4*i+3]];
    }
    //check if there are any remaining characters and pad with 0
    char leftover=sequence_length%4;
    if (!leftover){
        return last_element;
    }
    switch (leftover){
        case 1:
            code[last_element]=seq2code[(unsigned char)string[4*last_element]]<<6;
            break;
        case 2:
            code[last_element]=seq2code[(unsigned char)string[4*last_element]]<<6 | seq2code[(unsigned char)string[4*last_element+1]]<<4;
            break;
        case 3:
            code[last_element]=seq2code[(unsigned char)string[4*last_element]]<<6 | seq2code[(unsigned char)string[4*last_element+1]]<<4 | seq2code[(unsigned char)string[4*last_element+2]]<<2;
            break;
    }
    return last_element+1;
}
void string2all_codes(char *string, unsigned char codes[][LINE_LENGTH/2+1], int *lengths){
    //4 codes are returned for the string to capture all the possible frames
    char offset[3][LINE_LENGTH];
    const int seqlength = strlen(string);

    strcpy(offset[0], string+1);
    strcpy(offset[1], string+2);
    strcpy(offset[2], string+3);
    lengths[0]=string2code(string, seqlength, codes[0]);

    for (int i=0; i<3; i++){
        lengths[i+1]=string2code(offset[i], seqlength-1-i, codes[i+1]);
    }   
}

int check_barcode(char *barcode, unsigned char *code){
    for (int i=0; i<barcode_code_length; i++){
        string2code(barcode, barcode_length, code);
    }
    if (g_hash_table_lookup(whitelist_hash,code) == NULL){
        return 1;
    }
    return 0;
}
int find_matching_codes(unsigned char *sequence_code, unsigned char *feature_code, size_t feature_code_length, size_t feature_length){
    unsigned char a1[feature_code_length];
    unsigned char *b=feature_code;
    memcpy(a1, sequence_code, feature_code_length);
    unsigned char right_overhang=feature_length % 4;

    if (right_overhang){
        a1[feature_code_length-1] = a1[feature_code_length-1] & (0xff << (8-2*right_overhang));        
    }
    int length=feature_code_length;
    unsigned char *a=a1;
    if(length >= 8){
        uint64_t *a64=(uint64_t*) a;
        uint64_t *b64=(uint64_t*) b;
        for (int i=0; i<length/8; i++){
            if (a64[i] != b64[i]){
                return feature_code_length+1;
            }
        }
        if (length % 8 == 0){
            return 0;
        }
        a+=8*(length/8);
        b+=8*(length/8);
        length=length%8;
    }
    if(length >= 4){
        const uint32_t *a32=(uint32_t*) a;
        const uint32_t *b32=(uint32_t*) b;
        for (int i=0; i<length/4; i++){
            if (a32[i] != b32[i]){
                return feature_code_length+1;
            }
        }
        if (length % 4 == 0){
            return 0;
        }
        a+=4*(length/4);
        b+=4*(length/4);
        length=length%4;    
    }
    if (length >= 2){
        const uint16_t *a16=(uint16_t*) a;
        const uint16_t *b16=(uint16_t*) b;
        for (int i=0;i<length/2; i++){
            if (a16[i] != b16[i]){
                return feature_code_length+1;  
            }
        }
        if (length % 2 == 0){
            return 0;
        }
        a+=2*(length/2);
        b+=2*(length/2);
        length=length%2;
    }
    for (int i=0;i<length; i++){
        if (a[i] != b[i]){
            return length+1;
        }
    }
    return 0;
}
int fuzzy_matching_codes(unsigned char *sequence_code, unsigned char *feature_code, size_t feature_code_length, size_t feature_length, int maxHammingDistance){
    int mismatches=0;
    unsigned char a1[feature_code_length];
    unsigned char *b=feature_code;
    memcpy(a1, sequence_code, feature_code_length);
    unsigned char right_overhang=feature_length % 4;
    if (right_overhang){
        a1[feature_code_length-1] = a1[feature_code_length-1] & (0xff << (8-2*right_overhang));        
    }
    unsigned char *a=a1;
    int length=feature_code_length;
    if(length >= 8){
        const uint64_t *a64=(uint64_t*) a;
        const uint64_t *b64=(uint64_t*) b;
        for (int i=0; i<length/8; i++){
            if (a64[i] != b64[i]){
                if(mismatches==maxHammingDistance){
                    return feature_code_length+1;
                }
                uint64_t diff=a64[i]^b64[i];
                unsigned char *diff8=(unsigned char*)&diff;
                mismatches+=diff2Hamming[diff8[0]]+diff2Hamming[diff8[1]]+diff2Hamming[diff8[2]]+diff2Hamming[diff8[3]] + diff2Hamming[diff8[4]]+diff2Hamming[diff8[5]]+diff2Hamming[diff8[6]]+diff2Hamming[diff8[7]];
                if (mismatches > maxHammingDistance){
                    return feature_code_length+1;
                }
            }
        }
        if (length % 8 == 0){
            return mismatches;
        }
        a+=8*(length/8);
        b+=8*(length/8);
        length=length%8;
    }
    if (length >= 4 ){
        const int last=8*(length/8);
        const uint32_t *a32=(uint32_t*) a;
        const uint32_t *b32=(uint32_t*) b;
        for (int i=last/4; i<length/4; i++){
            if (a32[i] != b32[i]){
                if(mismatches==maxHammingDistance){
                    return feature_code_length+1;
                }
                uint32_t diff=a32[i]^b32[i];
                unsigned char *diff8=(unsigned char*)&diff;
                mismatches+=diff2Hamming[diff8[0]]+diff2Hamming[diff8[1]]+diff2Hamming[diff8[2]]+diff2Hamming[diff8[3]];
                if (mismatches > maxHammingDistance){
                    return feature_code_length+1;
                }
            }
        }
        if (length % 4 == 0){
            return mismatches;
        }
        a+=4*(length/4);
        b+=4*(length/4);
        length=length%4;
        
    }
    if (length >= 2){
        const uint16_t *a16=(uint16_t*) a;
        const uint16_t *b16=(uint16_t*) b;
        const int last=4*(length/4);
        for (int i=last/2;i<length/2; i++){
            if (a16[i] != b16[i]){
                if(mismatches==maxHammingDistance){
                    return feature_code_length+1;
                }
                uint16_t diff=a16[i]^b16[i];
                unsigned char *diff8=(unsigned char*)&diff;
                mismatches+=diff2Hamming[diff8[0]]+diff2Hamming[diff8[1]];  
                if (mismatches > maxHammingDistance){
                    return feature_code_length+1;
                }   
            }
        }
        if (length % 2 == 0){
            return mismatches;
        }
        a+=2*(length/2);
        b+=2*(length/2);
        length=length%2;
    }
    const int last=2*(length/2);
    for (int i=last;i<length; i++){
        if (a[i] != b[i]){
            if(mismatches==maxHammingDistance){
                return feature_code_length+1;
            }
            mismatches+=diff2Hamming[a[i]^b[i]];
            if (mismatches > maxHammingDistance){
                return feature_code_length+1;
            }
        }
    }   
    return mismatches;
}

int find_matches_in_sub_arrays(unsigned char *sequence_code, unsigned char *feature_code, size_t sequence_code_length, size_t feature_code_length, size_t feature_length, int maxHammingDistance, int *best_offset){
    int offset=0;
    int length= (sequence_code_length > feature_code_length) ? sequence_code_length : feature_code_length;
    int minHammingDistance=length+1;
    *best_offset=0;
    while (sequence_code_length - offset >= feature_code_length){
        if (!maxHammingDistance){
            if (!find_matching_codes(sequence_code+offset, feature_code, feature_code_length, feature_length)){
                *best_offset=offset;
                return 0;
            }
        }
        else{
            int hammingDistance=fuzzy_matching_codes(sequence_code+offset, feature_code, feature_code_length, feature_length, maxHammingDistance);
            if (hammingDistance < minHammingDistance){
                minHammingDistance=hammingDistance;
                *best_offset=offset;
            }
            //if we have found a perfect match return immediately
            if  (!minHammingDistance){
                *best_offset=offset;
                return 0;
            }   
        }
        offset++;
    }
    return minHammingDistance;
}

int checkSequenceAndCorrectForN(char *line, char *corrected_lines[], char *buffer,int sequence_length, int maxN){
    int nCount=0;
    int indices[maxN];
    corrected_lines[0]=line;
    for (int i=0; i<sequence_length; i++){
        if (line[i] == 'N'){
            indices[nCount++]=i;
            if (nCount >=  maxN) return 0;
            
        }
        else if (line[i] != 'A' && line[i] != 'C' && line[i] != 'G' && line[i] != 'T'){
            return 0;
        }
    }
    if (!nCount){
        return 1;
    }
    char line_copy[LINE_LENGTH];
    strcpy(line_copy, line);
    return generate_sequences(line_copy, sequence_length, indices, &buffer, 0, nCount, corrected_lines, 0);
}

char* printToBuffer(char *string, int sequence_length, char *buffer){
    for (int i=0; i<sequence_length; i++){
        buffer[i]=string[i];
    }
    buffer[sequence_length]='\0';
    return buffer + sequence_length+1;
}

int generate_sequences(char* string, int string_length, int *indices, char **output, int index, int number_of_indices, char *output_indices[], 
int number_of_outputs) {
    static char bases[] = "ACGT";
    //base case when index is equal to number of indices
    if (index == number_of_indices) {
        printToBuffer(string, string_length, *output);
        output_indices[number_of_outputs++]=*output;
        *output += string_length + 1;
        return number_of_outputs;
    }
    //when the index is less than the number of indices
    for (int i = 0; i < 4; i++) {
        string[indices[index]] = bases[i];  
        number_of_outputs=generate_sequences(string, string_length, indices, output, index+1, number_of_indices, output_indices, number_of_outputs);
    }
    return number_of_outputs;
}
int simple_search(feature_arrays *features, char *line){  
    for (int i=0; i<features->number_of_features; i++){
        char *query=line;
        char *feature=features->feature_sequences[i];
        //DEBUG_PRINT( "Comparing %s %s\n", query, feature);
        while (*query == *feature){
            query++;
            feature++;
        }
        if (!*feature){
            return i+1;
        }
    }
    return 0;
}
int simple_hamming_search(feature_arrays *features, char *line, int maxHammingDistance, int *hamming_distance){  
    int ambiguous=0;
    int bestFeature=0;
    int bestHammingDistance=maxHammingDistance;
    for (int i=0; i<features->number_of_features; i++){
        char *query=line;
        char *feature=features->feature_sequences[i];
        int hammingDistance=0;
        
        //DEBUG_PRINT( "Comparing %s %s\n", query, feature);
        while (*query && *feature){
            if (*query != *feature){
                hammingDistance++;
                if (hammingDistance > maxHammingDistance){
                    break;
                }
            }
            query++;
            feature++;
        }
        if (!*feature){
            if (hammingDistance == bestHammingDistance){
                if(bestFeature){
                    ambiguous=1;
                }
                else{
                    bestFeature=i+1;
                    bestHammingDistance=hammingDistance;
                    ambiguous=0;
                }
            }
            else if (hammingDistance < bestHammingDistance){
               bestFeature=i+1;
               bestHammingDistance=hammingDistance;
               ambiguous=0;
            }
        }

    }
    *hamming_distance=bestHammingDistance;
    if (bestHammingDistance <= maxHammingDistance && !ambiguous){
        return bestFeature;
    }
    return 0;
}
int find_feature_match_single(feature_arrays *features, char *lineR2, int maxHammingDistance,int *bestScore, char **matching_sequence, uint16_t *match_position){
    // convert lineR2 to 4 codes
    // do a quick check to see if there is a perfect match the constant feature
    //int best_feature=simple_search_code(features, lineR2, constant_offset);
    unsigned char codes[4][LINE_LENGTH/2+1];
    int code_lengths[4];
    string2all_codes(lineR2, codes, code_lengths);
    int best_feature=0;
    int bestHammingDistance=maxHammingDistance; 
    int best_sequence_offset=0;
    int ambiguous=0;
    for (int i=0; i<4; i++){
        for (int j=0; j<features->number_of_features; j++){
            int code_offset=0;
            int hammingDistance=find_matches_in_sub_arrays(codes[i], features->feature_codes[j], code_lengths[i], features->feature_code_lengths[j],features->feature_lengths[j], maxHammingDistance, &code_offset);
            if (!hammingDistance){
                *bestScore=0;
                best_feature=j+1;
                *matching_sequence=lineR2+i+code_offset*4;
                *match_position = i+code_offset*4;
                return best_feature;
            }
            if (hammingDistance < bestHammingDistance){
                best_feature=j+1;
                bestHammingDistance=hammingDistance;
                best_sequence_offset=i+4*code_offset;
                ambiguous=0;
            }
            else if (hammingDistance == bestHammingDistance){
                ambiguous=1;
            }
        }
    }
    if (bestHammingDistance <= maxHammingDistance){
        *bestScore=bestHammingDistance;
        *matching_sequence=lineR2+best_sequence_offset;
        if (!ambiguous){
            return best_feature;
        }
        return 0;
    }
    return 0;
}
int find_feature_match_parallel(feature_arrays *features, char *lineR2, int maxHammingDistance, int nThreads, int *bestScore, char **matching_sequence, uint16_t *match_position){
    // convert lineR2 to 4 codes
    // do a quick check to see if there is a perfect match the constant feature
    //int best_feature=simple_search_code(features, lineR2, constant_offset);
    nThreads=(nThreads > 4) ? 4 : nThreads;
    nThreads=(nThreads < 1) ? 1 : nThreads;
    if (nThreads==1){
        return find_feature_match_single(features, lineR2, maxHammingDistance, bestScore, matching_sequence, match_position);
    }
    unsigned char codes[4][LINE_LENGTH/2+1];
    int code_lengths[4];
    string2all_codes(lineR2, codes, code_lengths);

    int best_feature=0;
    int bestFeatureDistance=maxHammingDistance;   
    int bestHammingDistances[4]={maxHammingDistance,maxHammingDistance,maxHammingDistance,maxHammingDistance};
    int best_code_offsets[4]={0,0,0,0};
    int ambiguous[4]={0,0,0,0};
    int best_match[4]={0,0,0,0};
    int best_matching_indices[4]={0,0,0,0};

    int exact_match_found=0;
    #pragma omp parallel for num_threads(nThreads)
    for (int i=0; i<4; i++){
        //get the thread number
        const int thread_num=omp_get_thread_num();
        for (int j=0; j<features->number_of_features && !exact_match_found; j++){
            int code_offset=0;
            int hammingDistance=find_matches_in_sub_arrays(codes[i], features->feature_codes[j], code_lengths[i], features->feature_code_lengths[j],features->feature_lengths[j], maxHammingDistance, &code_offset);
            if (!hammingDistance){
                *bestScore=0;
                exact_match_found=1;
                best_feature=j+1;
                *matching_sequence=lineR2+i+code_offset*4;
                *match_position = i+code_offset*4;
                best_matching_indices[thread_num]=i;
                break;
            }
            if (hammingDistance < bestHammingDistances[thread_num]){
                best_match[thread_num]=j+1;
                bestHammingDistances[thread_num]=hammingDistance;
                best_code_offsets[thread_num]=code_offset;
                ambiguous[thread_num]=0;
                best_matching_indices[thread_num]=i;
            }
            else if (hammingDistance == bestHammingDistances[thread_num]){
                if (best_match[thread_num]){
                    ambiguous[thread_num]=1;
                }
                else{
                    ambiguous[thread_num]=0;
                    best_match[thread_num]=j+1;
                    best_code_offsets[thread_num]=code_offset;
                    bestHammingDistances[thread_num]=hammingDistance;
                    best_matching_indices[thread_num]=i;
                }
            }
        }
    }
    if (exact_match_found){//matching_sequence is already set
        return best_feature;
    }
    //find the best hamming distance - if there is a perfect match return immediately
    char multiAmbiguous=0;
    int best_code_offset=0;
    int best_i=0;
    for (int i=0; i<nThreads; i++){
        if (bestHammingDistances[i] < bestFeatureDistance){
            bestFeatureDistance=bestHammingDistances[i];
            best_feature=best_match[i];
            best_code_offset=best_code_offsets[i];
            multiAmbiguous=ambiguous[i]; //reset ambiguous flag
            best_i=best_matching_indices[i];
        }
        else if (bestHammingDistances[i] == bestFeatureDistance){
            if (best_feature ){
                multiAmbiguous=1;
            }
            else{
                best_feature=best_match[i];
                bestFeatureDistance=bestHammingDistances[i];
                best_code_offset=best_code_offsets[i];
                best_i=best_matching_indices[i];
                multiAmbiguous=ambiguous[i];
            }
        }
    }
    if (bestFeatureDistance <= maxHammingDistance){
        *bestScore=bestFeatureDistance;
        *matching_sequence=lineR2+best_i+best_code_offset*4;
        *match_position=best_i+best_code_offset*4;
        if (!multiAmbiguous){
            return best_feature;
        }
        return 0;
    }
    return 0;
}
void update_feature_counts(char *barcodeString, char *umi, uint32_t feature_index, data_structures *hashes,  memory_pool_collection *pools){ 
    unsigned char code[barcode_code_length];
    string2code(barcodeString, barcode_length, code);
    update_feature_counts_from_code(code, umi, feature_index, hashes, pools);
}
// In src/assignBarcodes.c

void update_feature_counts_from_code(unsigned char *code, char *umi, uint32_t feature_index, data_structures *hashes, memory_pool_collection *pools){
    update_umi_counts(code, umi, feature_index, hashes, pools);
    feature_counts *s = g_hash_table_lookup(hashes->filtered_hash, code);

    if (s == NULL) {
        // Logic for a NEW barcode
        s = (feature_counts*) allocate_memory_from_pool(pools->feature_counts_pool);
        if (s == NULL) return; // Error check

        memcpy(s->sequence_code, code, barcode_code_length);

        // NEW: Initialize the inner hash table for this barcode's counts.
        // We don't provide destroy functions because the keys/values are not allocated pointers.
        s->counts = g_hash_table_new(g_direct_hash, g_direct_equal);
        if (s->counts == NULL) {
            fprintf(stderr, "Fatal: GHashTable creation for counts failed.\n");
            exit(EXIT_FAILURE);
        }

        g_hash_table_insert(hashes->filtered_hash, s->sequence_code, s);

        // For a new entry, the count for this feature is 1.
        g_hash_table_insert(s->counts, GUINT_TO_POINTER(feature_index), GUINT_TO_POINTER(1));
        // The total count (at index 0) is also 1.
        g_hash_table_insert(s->counts, GUINT_TO_POINTER(0), GUINT_TO_POINTER(1));

    } else {
        // Logic for an EXISTING barcode
        // 1. Increment the specific feature's count
        uintptr_t current_count = GPOINTER_TO_UINT(g_hash_table_lookup(s->counts, GUINT_TO_POINTER(feature_index)));
        current_count++;
        g_hash_table_replace(s->counts, GUINT_TO_POINTER(feature_index), GUINT_TO_POINTER(current_count));

        // 2. Increment the total count stored at key '0'
        uintptr_t total = GPOINTER_TO_UINT(g_hash_table_lookup(s->counts, GUINT_TO_POINTER(0)));
        total++;
        g_hash_table_replace(s->counts, GUINT_TO_POINTER(0), GUINT_TO_POINTER(total));
    }
}
void update_umi_counts(unsigned char *code, char *umi,  uint32_t feature_index,data_structures *hashes, memory_pool_collection *pools){
    unsigned char code8[8];
    memset(code8,0,8);
    memcpy(code8,code,barcode_code_length);
    string2code(umi, umi_length, code8+barcode_code_length);

    feature_umi_counts *s = g_hash_table_lookup(hashes->sequence_umi_hash, code8);
    
    if (s == NULL) {
        // Logic for a NEW barcode-UMI combination
        s = (feature_umi_counts*) allocate_memory_from_pool(pools->feature_umi_counts_pool);
        if (s == NULL) return; // Error check

        memcpy(s->sequence_umi_code, code8, 8);
        
        // Create the new inner hash table for this UMI's feature counts
        s->counts = g_hash_table_new(g_direct_hash, g_direct_equal);
        
        // Set initial count for this feature to 1
        g_hash_table_insert(s->counts, GUINT_TO_POINTER(feature_index), GUINT_TO_POINTER(1));
        // Mark as unvisited for the connected components algorithm by setting key '0' to 0.
        g_hash_table_insert(s->counts, GUINT_TO_POINTER(0), GUINT_TO_POINTER(0));

        g_hash_table_insert(hashes->sequence_umi_hash, s->sequence_umi_code, s);
    } else {
        // Logic for an EXISTING barcode-UMI combination
        uintptr_t current_count = GPOINTER_TO_UINT(g_hash_table_lookup(s->counts, GUINT_TO_POINTER(feature_index)));
        current_count++;
        g_hash_table_replace(s->counts, GUINT_TO_POINTER(feature_index), GUINT_TO_POINTER(current_count));
    }
}
char check_neighbor(uint64_t code64,uint32_t *counts, data_structures *hashes){
    feature_umi_counts *result = g_hash_table_lookup(hashes->sequence_umi_hash, &code64);

    if (result && result->counts) {
        // Check if this node has been visited using key '0'
        uintptr_t visited_flag = GPOINTER_TO_UINT(g_hash_table_lookup(result->counts, GUINT_TO_POINTER(0)));
        
        if (visited_flag == 0) { // If not visited
            // --- NEW: Iterate through the hash table to sum counts ---
            GHashTableIter iter;
            gpointer key, value;
            g_hash_table_iter_init(&iter, result->counts);
            while (g_hash_table_iter_next(&iter, &key, &value)) {
                uint32_t feature_index = GPOINTER_TO_UINT(key);
                if (feature_index > 0) { // Don't add the visited flag to the temp counts array
                    counts[feature_index] += GPOINTER_TO_UINT(value);
                }
            }

            // Mark as visited by replacing the value at key '0' with 1
            g_hash_table_replace(result->counts, GUINT_TO_POINTER(0), GUINT_TO_POINTER(1));
            return 1; // Indicate success
        }
    }
    return 0; // Not found or already visited
}
int find_neighbors(uint64_t key64, uint64_t *neighbors, uint32_t *counts, data_structures *hashes){
    int neighbor_count=0;
    uint64_t code64=key64;
    unsigned char *code8= (unsigned char*) &code64;
    unsigned char *code=code8+barcode_code_length;
    for (int i=0; i<umi_code_length; i++){       
        const unsigned char mask=code[i];
        for (unsigned char j=1; j<4; j++){
            code[i]=mask ^ j;
            if (check_neighbor(code64,counts,hashes)) neighbors[neighbor_count++]=code64;
            code[i]=mask ^ (j << 2);
            if (check_neighbor(code64,counts,hashes)) neighbors[neighbor_count++]=code64;
            code[i]=mask ^ (j << 4);
            if (check_neighbor(code64,counts,hashes)) neighbors[neighbor_count++]=code64;
            code[i]=mask;
        }
    }
    return neighbor_count;
}
void find_deduped_counts(data_structures *hashes, GHashTable* barcode_to_deduped_counts, uint16_t stringency, uint16_t min_counts){
    GHashTableIter iter;
    gpointer lookup_key, result;
    g_hash_table_iter_init(&iter, hashes->sequence_umi_hash);
    uint32_t clique_counts[number_of_features+1];

    while (g_hash_table_iter_next(&iter, &lookup_key, &result)) {
        feature_umi_counts *umi_counts = (feature_umi_counts*) result;
        if (GPOINTER_TO_UINT(g_hash_table_lookup(umi_counts->counts, GUINT_TO_POINTER(0)))) {
            continue; // Already processed if the flag is 1
        }
        memset(clique_counts, 0, sizeof(clique_counts));
        
        find_connected_component(lookup_key, clique_counts, hashes);

        feature_counts *s = g_hash_table_lookup(hashes->filtered_hash, umi_counts->sequence_umi_code);
        if (s) {
            // Get or create the inner hash table for this barcode's deduped counts
            GHashTable* temp_deduped_hash = g_hash_table_lookup(barcode_to_deduped_counts, s->sequence_code);
            if (temp_deduped_hash == NULL) {
                temp_deduped_hash = g_hash_table_new(g_direct_hash, g_direct_equal);
                g_hash_table_insert(barcode_to_deduped_counts, s->sequence_code, temp_deduped_hash);
            }
            
            add_deduped_count(temp_deduped_hash, clique_counts, stringency, min_counts);
        }
    }
}

void find_connected_component(gpointer start_key, uint32_t *counts, data_structures *hashes){
    uint64_t neighbors[umi_length*3];
    clear_queue(hashes->neighbors_queue);
    enqueue(hashes->neighbors_queue,*(uint64_t*)start_key);
    while (!is_empty(hashes->neighbors_queue)) {
        uint64_t code = dequeue(hashes->neighbors_queue);
        check_neighbor(code,counts,hashes);
        int neighbor_count=find_neighbors(code,neighbors, counts, hashes);
        //if (neighbor_count) DEBUG_PRINT( "Found %d neighbors\n",neighbor_count);
        for (int i = 0; i < neighbor_count; i++) {
            enqueue(hashes->neighbors_queue,neighbors[i]);
        }
    }
}
/**
 * @brief Calculates the winning feature based on UMI clique counts and stringency,
 * and updates a temporary hash table with the deduplicated count.
 *
 * @param temp_deduped_hash A temporary hash table mapping a feature_index to its deduped count.
 * @param clique_counts An array of raw UMI counts for the connected component.
 * @param stringency The user-defined stringency for deduplication.
 * @param min_counts The minimum count threshold.
 */
 void add_deduped_count(GHashTable* temp_deduped_hash, uint32_t *clique_counts, uint16_t stringency, uint16_t min_counts) {
    uint32_t winning_feature_index = 0;

    // --- STRINGENCY LOGIC (Restored from original code) ---
    if (!stringency) {
        // RNA-seq strategy: any feature with enough counts gets a single deduped count.
        // This is a special case as it can increment multiple features.
        for (int i = 1; i < number_of_features + 1; i++) {
            if (clique_counts[i] > min_counts) {
                uintptr_t count = GPOINTER_TO_UINT(g_hash_table_lookup(temp_deduped_hash, GUINT_TO_POINTER(i)));
                count++;
                g_hash_table_replace(temp_deduped_hash, GUINT_TO_POINTER(i), GUINT_TO_POINTER(count));
            }
        }
        return; // Return early as the logic is different
    }

    if (stringency >= 1000) {
        // Highest stringency: only one feature can have counts.
        uint32_t feature_index = 0;
        uint32_t total_counts = 0;
        for (int i = 1; i < number_of_features + 1; i++) {
            if (clique_counts[i]) {
                if (feature_index) { // If a feature has already been found, there's more than one.
                    winning_feature_index = 0; // Invalidate winner
                    break;
                }
                feature_index = i;
                total_counts = clique_counts[i];
            }
        }
        if (feature_index && total_counts > min_counts) {
            winning_feature_index = feature_index;
        }
    } else {
        // Mid-stringency: find the feature with the highest, unique count.
        uint32_t feature_index = 0;
        uint32_t max_counts = 0;
        uint32_t total_counts = 0;
        unsigned char unique = 0;
        for (int i = 1; i < number_of_features + 1; i++) {
            if (clique_counts[i]) {
                total_counts += clique_counts[i];
                if (clique_counts[i] > max_counts) {
                    max_counts = clique_counts[i];
                    feature_index = i;
                    unique = 1;
                } else if (clique_counts[i] == max_counts) {
                    unique = 0;
                }
            }
        }
        if (unique && (double)max_counts > total_counts * stringency / 1000.0 && total_counts > min_counts) {
            winning_feature_index = feature_index;
        }
    }

    // --- If a single winner was found, increment its count in the temporary hash table ---
    if (winning_feature_index > 0) {
        uintptr_t current_deduped_count = GPOINTER_TO_UINT(
            g_hash_table_lookup(temp_deduped_hash, GUINT_TO_POINTER(winning_feature_index))
        );
        current_deduped_count++;
        g_hash_table_replace(temp_deduped_hash, 
                             GUINT_TO_POINTER(winning_feature_index), 
                             GUINT_TO_POINTER(current_deduped_count));
    }
}

void code2string(unsigned char *code, char *string, int length){
    for (int i=0; i<length; i++){
        string[4*i]=code2seq[code[i]][0];
        string[4*i+1]=code2seq[code[i]][1];
        string[4*i+2]=code2seq[code[i]][2];
        string[4*i+3]=code2seq[code[i]][3];
    }
    string[4*length]='\0';
}

void printFeatureCounts(feature_arrays *features, int *deduped_counts, int *barcoded_counts,int **coexpression_counts, int **coexpression_histograms, char *directory, data_structures *hashes, statistics *stats, uint16_t stringency, uint16_t min_counts){
    int total_deduped_counts = 0;
    int total_raw_counts = 0;
    char barcodes_file[FILENAME_LENGTH];
    char stats_file[FILENAME_LENGTH];

    memset(deduped_counts, 0, features->number_of_features * sizeof(int));
    memset(barcoded_counts, 0, features->number_of_features * sizeof(int));

    mkdir_p(directory);
    sprintf(barcodes_file, "%s/barcodes.txt", directory);
    sprintf(stats_file, "%s/stats.txt", directory);
    FILE *barcodesfp = fopen(barcodes_file, "w");
    if (barcodesfp == NULL) { /* ... error handling ... */ }

    char matrix_file[LINE_LENGTH];
    sprintf(matrix_file, "%s/features_matrix.mtx", directory);
    FILE *matrixfp = fopen(matrix_file, "w");
    if (matrixfp == NULL) { /* ... error handling ... */ }

    // --- Step 1: Create and populate the temporary hash for deduped counts ---
    //
    // CRITICAL FIX: Use the SAME hash and equal functions as filtered_hash
    // to ensure keys (the sequence_code pointers) are handled consistently.
    // Also, tell GLib how to destroy the inner hash tables that are used as values.
    //
    GHashTable* barcode_to_deduped_hash = g_hash_table_new_full(hash_int32, equal_int32, NULL, (GDestroyNotify)g_hash_table_destroy);

    find_deduped_counts(hashes, barcode_to_deduped_hash, stringency, min_counts);

    // --- Step 2: Calculate stats and write the Matrix Market header ---
    size_t number_of_features_seen = 0;
    size_t number_of_barcode_entries = g_hash_table_size(barcode_to_deduped_hash);
    GHashTableIter iter;
    gpointer key, value;
    g_hash_table_iter_init(&iter, barcode_to_deduped_hash);
    while (g_hash_table_iter_next(&iter, &key, &value)) {
        number_of_features_seen += g_hash_table_size((GHashTable*)value);
    }
    
    fprintf(matrixfp, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(matrixfp, "%%metadata_json: {\"software_version\": \"assignBarcodes-0.1\", \"format_version\": 1}\n");
    fprintf(matrixfp, "%d %ld %ld\n", features->number_of_features, number_of_barcode_entries, number_of_features_seen);
    
    
    
    // --- NEW (Restored Logic): Populate Co-expression Matrices ---
    // This loop must run before we start writing the output file.
    if (coexpression_counts != NULL && coexpression_histograms != NULL) {
        g_hash_table_iter_init(&iter, barcode_to_deduped_hash);
        while (g_hash_table_iter_next(&iter, &key, &value)) {
            GHashTable* deduped_hash = (GHashTable*)value;
            int number_of_coexpressors = 0;
            int coCounts[features->number_of_features + 1];
            int coExpressorsIndices[features->number_of_features];
            
            GHashTableIter inner_iter;
            gpointer inner_key, inner_value;
            g_hash_table_iter_init(&inner_iter, deduped_hash);
            while(g_hash_table_iter_next(&inner_iter, &inner_key, &inner_value)) {
                int feature_index = GPOINTER_TO_INT(inner_key);
                int count = GPOINTER_TO_INT(inner_value);
                if (count > 0) {
                    coCounts[number_of_coexpressors] = count;
                    coExpressorsIndices[number_of_coexpressors++] = feature_index;
                }
            }
            
            if (number_of_coexpressors > 1) {
                for (int i = 0; i < number_of_coexpressors; i++) {
                    for (int j = 0; j < number_of_coexpressors; j++) {
                        coexpression_counts[coExpressorsIndices[i]][coExpressorsIndices[j]] += coCounts[i];
                    }
                    coexpression_histograms[coExpressorsIndices[i]][number_of_coexpressors]++;
                }
            } else if (number_of_coexpressors) {
                coexpression_counts[coExpressorsIndices[0]][0] += coCounts[0];
                coexpression_counts[0][coExpressorsIndices[0]] += coCounts[0];
                coexpression_histograms[coExpressorsIndices[0]][1]++;
            }
        }
    }
    // --- Step 3: Iterate through all barcodes to write files and calculate final stats ---
    int line_no = 1;
    g_hash_table_iter_init(&iter, hashes->filtered_hash);
    while (g_hash_table_iter_next(&iter, &key, &value)) {
        feature_counts *entry = (feature_counts*)value;
        total_raw_counts += GPOINTER_TO_UINT(g_hash_table_lookup(entry->counts, GUINT_TO_POINTER(0)));
        // --- Populate total_barcoded_counts (raw counts) ---
        GHashTableIter raw_iter;
        gpointer raw_key, raw_value;
        g_hash_table_iter_init(&raw_iter, entry->counts);
        while (g_hash_table_iter_next(&raw_iter, &raw_key, &raw_value)) {
            int feature_index = GPOINTER_TO_INT(raw_key);
            if (feature_index > 0) { // Skip the total stored at index 0
                barcoded_counts[feature_index - 1] += GPOINTER_TO_INT(raw_value);
            }
        }

        // Check if this barcode has any deduplicated counts
        GHashTable* deduped_hash = g_hash_table_lookup(barcode_to_deduped_hash, entry->sequence_code);
        if (deduped_hash && g_hash_table_size(deduped_hash) > 0) {
            
            // Write barcode string to barcodes.txt
            char barcode[barcode_length + 1];
            code2string(entry->sequence_code, barcode, barcode_code_length);
            fprintf(barcodesfp, "%s\n", barcode);
                    // --- Populate total_barcoded_counts (raw counts) ---
        GHashTableIter raw_iter;
        gpointer raw_key, raw_value;
        g_hash_table_iter_init(&raw_iter, entry->counts);
        while (g_hash_table_iter_next(&raw_iter, &raw_key, &raw_value)) {
            int feature_index = GPOINTER_TO_INT(raw_key);
            if (feature_index > 0) { // Skip the total stored at index 0
                barcoded_counts[feature_index - 1] += GPOINTER_TO_INT(raw_value);
            }
        }
            // Write its matrix entries
            GHashTableIter dedup_iter;
            gpointer dedup_key, dedup_value;
            g_hash_table_iter_init(&dedup_iter, deduped_hash);
            while (g_hash_table_iter_next(&dedup_iter, &dedup_key, &dedup_value)) {
                int deduped_count = GPOINTER_TO_INT(dedup_value);
                fprintf(matrixfp, "%d %d %d\n", GPOINTER_TO_INT(dedup_key), line_no, deduped_count);
                total_deduped_counts += deduped_count;
            }
            line_no++;
        }
        
        // Clean up the raw counts hash table for this entry
        g_hash_table_destroy(entry->counts);
    }
    
    // --- Step 4: Final Cleanup ---
    // This single call will now correctly destroy the outer hash table
    // AND trigger g_hash_table_destroy on each of the inner hash tables.
    g_hash_table_destroy(barcode_to_deduped_hash);
    
    fclose(barcodesfp);
    fclose(matrixfp);
    fprintf(stderr,"closing matrix file\n");
    fprintf(stderr,"writing stats file\n");
    FILE *statsfp = fopen(stats_file, "w");
    fprintf (stderr, "Total feature counts %d\n", total_raw_counts);
    fprintf (stderr, "Total deduped feature counts %d\n", total_deduped_counts);
    fprintf (stderr, "Total unique barcode UMIs %d\n", g_hash_table_size(hashes->sequence_umi_hash));
    fprintf (stderr, "Total whitelisted barcodes %d\n", g_hash_table_size(hashes->filtered_hash));
    fprintf (stderr,"Total feature counts %d total_unmatched_reads %ld\n", total_raw_counts, stats->total_unmatched_features);
    fprintf (stderr, "Percentage reads assigned to barcode %.4f\n", 100.0*(total_raw_counts/(double) (total_raw_counts+stats->total_unmatched_features)));
    fprintf (statsfp, "Total feature counts %d\n", total_raw_counts);
    fprintf (statsfp, "Total deduped feature counts %d\n", total_deduped_counts);
    fprintf (statsfp, "Total unique barcode UMIs %d\n", g_hash_table_size(hashes->sequence_umi_hash));
    fprintf (statsfp, "Total whitelisted barcodes %d\n", g_hash_table_size(hashes->filtered_hash));
    fprintf (statsfp,"Total_unmatched_reads %ld\n", stats->total_unmatched_features);
    fprintf (statsfp, "Percentage reads assigned to barcode %.4f\n", 100.0*(total_raw_counts/(double) (total_raw_counts+stats->total_unmatched_features)));
    fclose(statsfp);

}
int find_closest_barcodes(unsigned char* code,unsigned char *corrected_codes, unsigned char *indices){
    int number_of_variants=0;
    for (int i=0; i<barcode_length; i++){
        unsigned char corrected_bases[3];
        //check whitelist for variant match
        int nmatches=find_variant_match(code, i, corrected_bases);
        if(nmatches){
            if (nmatches+number_of_variants > max_barcode_mismatches){
                return 0;
            }
            for (int j=0; j<nmatches && number_of_variants+j < max_barcode_mismatches; j++){
                memcpy(corrected_codes+(number_of_variants+j)*barcode_code_length, code, barcode_code_length);
                const int code_index = i / 4;
                corrected_codes[(number_of_variants+j)*barcode_code_length + code_index] = corrected_bases[j];
                indices[number_of_variants+j]=i;
            }
            number_of_variants+=nmatches;
        }
        if (number_of_variants > max_barcode_mismatches){
            return 0;
        }
    }
    return number_of_variants;
}
int find_variant_match(unsigned char *code, int sequence_index, unsigned char *corrected_bases ){
    const unsigned char index= sequence_index / 4;
    const unsigned char shift= 6-2*(sequence_index % 4);
    unsigned char mod_code[barcode_code_length];
    memcpy(mod_code, code, barcode_code_length);
    int number_of_variants=0; 
    for (unsigned char i=1; i<4; i++){
        unsigned char mod = code[index] ^ (i << shift );
        mod_code[index]=mod;
        if (g_hash_table_lookup(whitelist_hash, mod_code) != NULL){
            corrected_bases[number_of_variants++]=mod;
        }
    }
    return number_of_variants;
}
void process_pending_barcodes( data_structures *hashes, memory_pool_collection *pools, statistics *stats, double min_posterior){
    unmatched_barcodes_features_block *current_entry_block=stats->unmatched_list.first_entry;
    while (current_entry_block != NULL){
        unsigned char *retcode=find_best_posterior_match(current_entry_block, number_of_features, min_posterior,stats, hashes);
        if (retcode != 0){
            unmatched_barcodes_features current_entry;
            read_unmatched_features_block(current_entry_block, &current_entry);
            char umi[umi_length+1];
            code2string(current_entry.umi,umi, umi_code_length);
            update_feature_counts_from_code(retcode, umi, current_entry.feature_index, hashes,pools);
        }
        current_entry_block=current_entry_block->next;
    }
}
unsigned char* find_best_posterior_match (unmatched_barcodes_features_block *entry_block, int number_of_features, double min_posterior, statistics *stats, data_structures *hashes){
    unmatched_barcodes_features entry_struct;
    unmatched_barcodes_features *entry=&entry_struct;
    read_unmatched_features_block(entry_block, entry);
    if (!entry->feature_index){
        return 0;
    }
    double priors[max_barcode_mismatches+1];
    double evidence[max_barcode_mismatches+1];
    double total_evidence = 0.0;
    for (int i=0; i<entry->number_of_closest_barcodes; i++){
        //find counts for the barcode
        int total_counts=1;
        unsigned char* barcode =(entry->closest_barcodes) + i*barcode_code_length;
        feature_counts *s=g_hash_table_lookup(hashes->filtered_hash, barcode);
        if (s != NULL){
            total_counts += GPOINTER_TO_UINT(g_hash_table_lookup(s->counts, GUINT_TO_POINTER(0)));
        }
        priors[i]= pow(10,-0.1*(entry->Qscores[i]-33));
        evidence[i]=(double) total_counts * priors[i];
        total_evidence+=evidence[i];
    }
    for (int i = 0; i < entry->number_of_closest_barcodes; i++) {
        double posteriors = evidence[i] / total_evidence;
        if(posteriors > min_posterior){
            stats->pending_recovered++;
            return (entry->closest_barcodes) + i*barcode_code_length;
        }   
    }
    return 0;
}
int simpleCorrectFeature(char *line, feature_arrays *features, int maxN, int maxHammingDistance, int *hamming_distance){
    int best_feature=simple_search(features, line);
    *hamming_distance=0;
    if (!best_feature){
        const size_t length=strlen(line)-1;
        char buffer[(length+1) * (4 << ((maxN-1)*2))];
        char *corrected_seqs[ 4 << ((maxN-1)*2)];
        int nAlts=checkSequenceAndCorrectForN(line, corrected_seqs, buffer, length, maxN);
        if (nAlts == 1){
            return simple_hamming_search(features, corrected_seqs[0], maxHammingDistance,hamming_distance);
        }
        if (nAlts > 1){
            int bestHammingDistance=maxHammingDistance;
            int ambiguous=0;
            for (int i=0; i<nAlts; i++){
                int feature_index=simple_search(features, corrected_seqs[i]);
                if (feature_index){
                    *hamming_distance=0;
                    return feature_index;
                }
            }
            for (int i=0; i<nAlts; i++){
                int myHammingDistance;
                int feature_index=simple_hamming_search(features, corrected_seqs[i], maxHammingDistance,&myHammingDistance);
                if (feature_index && myHammingDistance < bestHammingDistance){
                    bestHammingDistance=myHammingDistance;
                    best_feature=feature_index;
                    ambiguous=0;
                }
                else if (feature_index && myHammingDistance == bestHammingDistance){
                    if(best_feature){
                        ambiguous=1;
                    }
                    else{
                        best_feature=feature_index;
                        bestHammingDistance=myHammingDistance;
                        ambiguous=0;
                    }
                }
            }
            if(best_feature && !ambiguous){
                *hamming_distance=bestHammingDistance;
                return best_feature;
            }
        }

    }
    return best_feature;
}
int checkAndCorrectFeature(char *line, feature_arrays *features,int maxHammingDistance, int nThreads, int *hamming_distance, char *matching_sequence, int maxN,char *ambiguous, uint16_t *match_position){
    const size_t length=strlen(line)-1;
    char buffer[(length+1) * (4 << ((maxN-1)*2))];
    char *corrected_seqs[ 4 << ((maxN-1)*2)];
    int hamming=0;
    //return ambiguous if the hamming distance is non-zero but the feature is zero
    //return ambiguous if there are too many Ns ie. nAlts is zero - distinguish this by setting hamming distance to maxHammingDistance+1
    int nAlts=checkSequenceAndCorrectForN(line, corrected_seqs, buffer, length, maxN);
    if (!nAlts){
        *hamming_distance=maxHammingDistance+1;
        *ambiguous=1;
        return 0;
    }
    if (nAlts == 1){
        char *myMatchingSequence;
        uint16_t myMatchPosition;
        //no Ns and the barcode is good in terms of ACGT 
        int feature_index=find_feature_match_parallel(features, line, maxHammingDistance,nThreads,&hamming,&myMatchingSequence, &myMatchPosition);
        *hamming_distance=hamming;
        if (feature_index){
            memcpy(matching_sequence, myMatchingSequence, features->feature_lengths[feature_index-1]);
            matching_sequence[features->feature_lengths[feature_index-1]]='\0';
            *match_position = myMatchPosition;
        }
        else if (hamming <= maxHammingDistance){
            *ambiguous=1;
        }    
        return feature_index;
    }
    else{
        //if there are Ns in the barcode then we need to check all the possible sequences
        //and return if there is a unique match in the whitelist
        int myAmbiguous=0;
        int bestfeature_index=0;
        int bestHammingDistance=maxHammingDistance;
        char *bestMatchingSequence=0; 
        uint16_t bestMatchPosition=0;
        //int bestAlt=0;

        for (int i=0; i<nAlts; i++){
            //DEBUG_PRINT( "Checking %s\n", corrected_seqs[i]);
            char *myMatchingSequence;
            uint16_t myMatchPosition;
            int feature_index=find_feature_match_parallel(features, corrected_seqs[i], maxHammingDistance,nThreads,&hamming,&myMatchingSequence, &myMatchPosition);
            //DEBUG_PRINT( "Feature index %d %d\n", feature_index,hamming);
            
            if (feature_index && !hamming){
                memcpy(matching_sequence, myMatchingSequence, features->feature_lengths[feature_index-1]);
                matching_sequence[features->feature_lengths[feature_index-1]]='\0';
                *hamming_distance=hamming;
                *match_position = myMatchPosition;
                return feature_index;
            }
            else if (hamming < bestHammingDistance){
                    bestHammingDistance=hamming;
                    bestfeature_index=feature_index;
                    bestMatchingSequence=myMatchingSequence;
                    myAmbiguous=0;
                    bestMatchPosition=myMatchPosition;
            }
            else if (hamming == bestHammingDistance){
                if( bestfeature_index){
                    myAmbiguous=1;
                }
            }
        }
        if (bestHammingDistance <= maxHammingDistance){
            *hamming_distance=bestHammingDistance;
            *ambiguous=myAmbiguous;
            if (!myAmbiguous && bestfeature_index){
                memcpy(matching_sequence, bestMatchingSequence, features->feature_lengths[bestfeature_index-1]);
                matching_sequence[features->feature_lengths[bestfeature_index-1]]='\0';
                *match_position = bestMatchPosition;
                return bestfeature_index;
            }
            return 0;
        }
        return 0;
    }
}
size_t barcode_code2number(unsigned char *code){
    unsigned char *key=(unsigned char*)g_hash_table_lookup(whitelist_hash, code);
    if (key == NULL){
        return 0;
    }
    return (key-whitelist)/barcode_code_length+1;
}
int checkAndCorrectBarcode(char **lines, int maxN, uint32_t feature_index, uint16_t match_position, data_structures *hashes, memory_pool_collection *pools, statistics *stats, int barcode_constant_offset){
    char *sequence=lines[0]+barcode_constant_offset;
    //DEBUG_PRINT( "Checking barcode %s\n", lines[1]);
    if (strlen(sequence) < barcode_length + umi_length){
        fprintf(stderr, "Error: Incomplete barcode %s\n", lines[0]);
        return 0;
    }
    if (!check_sequence(sequence + barcode_length, umi_length)){
        stats->nMismatches++;
        return 0;
    }
    //The return code should indicate whether we should calculate the barcode or not
    char buffer[(barcode_length + 1) * (4 << ((max_barcode_n-1)*2))];
    char *corrected_seqs[ 4 << ((max_barcode_n-1)*2)];
    char *candidateBarcode = sequence;
    memset (buffer, 0, (barcode_length + 1) * (4 << ((max_barcode_n-1)*2)));
 
    //DEBUG_PRINT( "Checking barcode %s\n", candidateBarcode);    
    int nAlts=checkSequenceAndCorrectForN(candidateBarcode, corrected_seqs, buffer, barcode_length, maxN);
    if (!nAlts){
        stats->nMismatches++;
        return 0;
    }
    unsigned char code[barcode_code_length];
    if (nAlts == 1){
        //no Ns and the barcode is good in terms of ACGT 
        //check if the barcode is in the filtered whitelist 
        string2code(candidateBarcode, barcode_length, code); 
        if (g_hash_table_lookup(whitelist_hash, code)){
            update_feature_counts_from_code(code, sequence+barcode_length, feature_index, hashes, pools);
            stats->valid++;            
            return 1;
        }
        //if the barcode is not in the whitelist then try to correct it 
        unsigned char indices[max_barcode_mismatches+1];
        unsigned char corrected_codes[(max_barcode_mismatches+1)*barcode_code_length];
        stats->nMismatches++;
        //try to correct it by finding closest barcodes
        int nMatches=find_closest_barcodes(code,corrected_codes, indices);
        if (!nMatches){
            return 0;
        }

        //if unique closest barcode found then correct
        if (nMatches==1){
            char corrected_barcode[barcode_length+1];   
            code2string(corrected_codes, corrected_barcode, barcode_code_length);
            memcpy(sequence, corrected_barcode, barcode_length);
            stats->recovered++;
            update_feature_counts_from_code(corrected_codes, sequence+barcode_length, feature_index, hashes, pools);
            stats->valid++;
            return 1;
        }
        else{
            //otherwise store it for later processing when we have counts to calculate the posterior probabilities
            if (nMatches > 1){
                unsigned char qscores[max_barcode_mismatches+1];
                for (int i=0; i<nMatches; i++){
                    qscores[i]=lines[1][indices[i]];
                }
                add_unmatched_barcode_store_feature(code, corrected_codes, sequence+barcode_length, qscores,feature_index,nMatches, match_position, pools,stats);
                stats->pending++;
            }
            return 0;
        }
        return 0;
    }
    else{
        stats->nMismatches++;
        //if there are Ns in the barcode then we need to check all the possible sequences
        //and return if there is a unique match in the whitelist
        int number_of_matches=0;
        unsigned char corrected_code[barcode_code_length];
        for (int i=0; i<nAlts; i++){
            unsigned char code[barcode_code_length];
            string2code(corrected_seqs[i], barcode_length, code);
            char corrected_barcode[barcode_length+1];   
            code2string(code, corrected_barcode, barcode_code_length);
            if (g_hash_table_lookup(whitelist_hash, code)){
                if (number_of_matches){
                    return 0;
                }
                else{
                    number_of_matches++;
                    memcpy(corrected_code, code, barcode_code_length);  
                }
            }
        }
        if (number_of_matches == 1){
            stats->recovered++;
            char corrected_barcode[barcode_length+1];   
            code2string(corrected_code, corrected_barcode, barcode_code_length);
            memcpy(sequence, corrected_barcode, barcode_length);
            update_feature_counts_from_code(corrected_code, sequence + barcode_length,feature_index, hashes,pools);
            stats->valid++;
            return 1;
        }
        return 0;
    }
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


void finalize_processing(feature_arrays *features, data_structures *hashes,  char *directory, memory_pool_collection *pools, statistics *stats, uint16_t stringency, uint16_t min_counts, double min_posterior){
    process_pending_barcodes(hashes, pools, stats,min_posterior);
    double elapsed_time = get_time_in_seconds() - stats->start_time;
    fprintf(stderr, "Finished processing %ld reads in %.2f seconds (%.1f thousand reads/second)\n", stats->number_of_reads, elapsed_time, stats->number_of_reads / (double)elapsed_time / 1000.0);

    //find size of whitelist_hash and non_whitelist_hash
    int total_feature_counts[features->number_of_features];
    int total_deduped_counts[features->number_of_features];
    int total_barcoded_counts[features->number_of_features];

    int **coExpression = NULL;
    int **coexpression_histograms = NULL;
    int *coExpressionStorage = NULL;
    int *histogramsStorage = NULL;
    //coexpression matrix will store the total counts for each feature (also when there are no others - in the 0 field) in cells with feature i starting with 1
    //for simplicity we also hava a zero line so that all the indices are 1 based
    //store the total counts in the zero line
    if (features->number_of_features >10000){
        printFeatureCounts(features, total_deduped_counts, total_barcoded_counts, coExpression, coexpression_histograms, directory, hashes, stats, stringency, min_counts);
        return;
    }
    coExpressionStorage = malloc((features->number_of_features + 1) * (features->number_of_features + 1) * sizeof(int));
    if (coExpressionStorage == NULL) {
        perror("Failed to allocate memory for coexpression matrix");
        exit(EXIT_FAILURE);
    }
    memset(coExpressionStorage, 0, (features->number_of_features + 1) * (features->number_of_features + 1) * sizeof(int));
    coExpression = malloc((features->number_of_features + 1) * sizeof(int *));
    if (coExpression == NULL) {
        perror("Failed to allocate memory for coexpression matrix");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < features->number_of_features + 1; i++) {
        coExpression[i] = coExpressionStorage + i * (features->number_of_features + 1);
    }
    histogramsStorage = malloc((features->number_of_features + 1) * (features->number_of_features + 1) * sizeof(int));
    if (histogramsStorage == NULL) {
        perror("Failed to allocate memory for coexpression histograms");
        exit(EXIT_FAILURE);
    }
    memset(histogramsStorage, 0, (features->number_of_features + 1) * (features->number_of_features + 1) * sizeof(int));
    coexpression_histograms = malloc((features->number_of_features + 1) * sizeof(int *));
    if (coexpression_histograms == NULL) {
        perror("Failed to allocate memory for coexpression histograms");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < features->number_of_features + 1; i++) {
        coexpression_histograms[i] = histogramsStorage + i * (features->number_of_features + 1);
    }
    //DEBUG_PRINT( "Number of reads matched to a feature %ld\n", valid);
    printFeatureCounts(features, total_deduped_counts, total_barcoded_counts, coExpression, coexpression_histograms, directory, hashes, stats, stringency, min_counts);
    //DEBUG_PRINT( "Percentage of reads matched to a feature %.2f\n", (valid / (double)number_of_reads * 100.0));

    DEBUG_PRINT( "Number of mismatched reads that are matched to a nearest barcode unambiguously  %ld\n", stats->recovered);
    DEBUG_PRINT( "Number of reads pending priors before matching to a barcode %ld\n", stats->pending);
    DEBUG_PRINT( "Number of pending reads that were successfully matched to a barcode %ld\n", stats->pending_recovered);
    
    print_feature_sequences(features, total_feature_counts, directory, hashes);
    char feature_stats_filename[FILENAME_LENGTH];
    sprintf(feature_stats_filename, "%s/features.txt", directory);
    FILE *feature_statsfp = fopen(feature_stats_filename, "w");
    if (feature_statsfp == NULL) {
        perror("Failed to open feature stats file");
        exit(EXIT_FAILURE);
    }
    fprintf(feature_statsfp, "Feature total_deduped_counts total_barcoded_counts total_feature_counts\n");

    char feature_coexpression_filename[FILENAME_LENGTH];
    sprintf(feature_coexpression_filename, "%s/feature_coexpression.txt", directory);
    FILE *feature_coexpression_fp = fopen(feature_coexpression_filename, "w");
    if (feature_coexpression_fp == NULL) {
        perror("Failed to open feature coexpression file");
        exit(EXIT_FAILURE);
    }
    int feature_printed[features->number_of_features];
    memset(feature_printed, 0, features->number_of_features * sizeof(int));
    int count = 0;
    for (int i = 0; i < features->number_of_features; i++) {
        if (feature_printed[i] == 0) {
            fprintf(feature_statsfp, "%s %d %d %d\n", features->feature_names[i], total_deduped_counts[count], total_barcoded_counts[count], total_feature_counts[count]);
            feature_printed[i] = 1;
            count++;
        }
    }
    for (int i = 0; i <= features->number_of_features; i++) {
        for (int j = 0; j < features->number_of_features + 1; j++) {
            fprintf(feature_coexpression_fp, "%d ", coExpression[i][j]);
        }
        fprintf(feature_coexpression_fp, "\n");
    }
    char feature_histograms_filename[FILENAME_LENGTH];
    sprintf(feature_histograms_filename, "%s/feature_histograms.txt", directory);
    FILE *feature_histograms_fp = fopen(feature_histograms_filename, "w");
    if (feature_histograms_fp == NULL) {
        perror("Failed to open feature histograms file");
        exit(EXIT_FAILURE);
    }
    //the 0 indices of the histograms to find the rightmost non-zero index
    for (int i = 1; i < features->number_of_features + 1; i++) {
        coexpression_histograms[i][0] = 0;
        for (int j = 1; j < features->number_of_features + 1; j++) {
            if (coexpression_histograms[i][j] > 0) {
                coexpression_histograms[i][0] = j;
            }
        }
    }
    for (int i = 1; i < features->number_of_features + 1; i++) {
        for (int j = 1; j < features->number_of_features + 1; j++) {
            fprintf(feature_histograms_fp, "%d ", coexpression_histograms[i][j]);
        }
        fprintf(feature_histograms_fp, "\n");
    }
    #ifndef NO_HEATMAP
    generate_heatmap(directory, features, coexpression_histograms);
    #endif
    free(coExpression);
    free(coExpressionStorage);
    free(coexpression_histograms);
    free(histogramsStorage);
    fclose(feature_statsfp);
    fclose(feature_coexpression_fp);
    fclose(feature_histograms_fp);
}

void read_fastq_chunks(gzFile barcode_fastqgz, gzFile forward_fastqgz, gzFile reverse_fastqgz, char *barcode_lines[], char *forward_lines[], char *reverse_lines[], char *buffer, char *done) {
    char *bufferptr = buffer;
    long long number_of_reads=0;
    gzFile gzfps[3] = {barcode_fastqgz, forward_fastqgz, reverse_fastqgz};  
    char **lines[3] = {barcode_lines, forward_lines, reverse_lines};
    for (int i=0; i<3; i++){
        gzFile gzfp=gzfps[i];
        if(!gzfp){
            continue;
        }
        fprintf(stderr, "Reading file %d %lld\n", i, number_of_reads);
        if (max_reads && number_of_reads >= max_reads){
            *done = 1;
            break;
        }
        for (int j = 0; j < 4; j++) {
            lines[i][j] = bufferptr;
            if (gzgets(gzfp, lines[i][j], LINE_LENGTH) != Z_NULL) {
            // Overwrite the newline character with null terminator
                const int length = strlen(lines[i][j]);
                if (lines[i][j][length - 1] != '\n') {
                    fprintf(stderr, "Error: Incomplete record in the FASTQ file\n");
                    exit(EXIT_FAILURE);
                }
                lines[i][j][length - 1] = '\0';
                if(i==2 && j==1){
                    if (j == 1) {
                        reverse_complement_in_place(bufferptr);
                    }
                    if (j == 3) {
                        reverse_in_place(bufferptr);
                    }
                }
                //DEBUG_PRINT( "Read %s\n", lines[i][j]);
                bufferptr += length;
            } 
            else {
                if (j == 0) {
                    *done = 1;
                    break;
                } else {
                 fprintf(stderr, "Error: Incomplete record in the FASTQ file\n");
                 exit(EXIT_FAILURE);
                }
            }
        }
        number_of_reads++;

    }    
}
void open_fastq_files(const char *barcode_fastq, const char *forward_fastq, const char *reverse_fastq, gzFile *barcode_fastqgz, gzFile *forward_fastqgz, gzFile *reverse_fastqgz) {
    *barcode_fastqgz = gzopen(barcode_fastq, "r");
    if (*barcode_fastqgz == NULL) {
        fprintf(stderr, "Error: Unable to open R1 FASTQ file %s\n", barcode_fastq);
        exit(EXIT_FAILURE);
    }
    if (forward_fastq != NULL) {
        *forward_fastqgz = gzopen(forward_fastq, "r");
        if (*forward_fastqgz == NULL) {
            fprintf(stderr, "Error: Unable to open R2 FASTQ file %s\n", forward_fastq);
            exit(EXIT_FAILURE);
        }
    }
    if (reverse_fastq != NULL) {
        *reverse_fastqgz = gzopen(reverse_fastq, "r");
        if (*reverse_fastqgz == NULL) {
            fprintf(stderr, "Error: Unable to open R3 FASTQ file %s\n", reverse_fastq);
            exit(EXIT_FAILURE);
        }
    }
}
fastq_reader* allocate_fastq_reader( char **filenames, int nfiles, int filetype, size_t read_size, size_t read_buffer_lines) {
    fastq_reader *reader = malloc(sizeof(fastq_reader));
    int total_filename_length = 0;
    for (int i = 0; i < nfiles; i++) {
        total_filename_length += strlen(filenames[i])+1; 
    }
    reader->filenames = malloc(nfiles * sizeof(char *));
    reader->concatenated_filenames=malloc(total_filename_length);

    char *this_filename = reader->concatenated_filenames;
    for (int i = 0; i < nfiles; i++) {
        fprintf (stderr, "Filename %s\n", filenames[i]);
        //do no increment the last time to avoid overun of the buffe
    }

    for (int i = 0; i < nfiles; i++) {
        strcpy(this_filename, filenames[i]);
        reader->filenames[i] = this_filename;
        this_filename += strlen(filenames[i]) + 1;
    }
    reader->gz_pointer = NULL;            // Initialize gz_pointer as needed
    reader->filetype = filetype;          // Set the filetype
    reader->read_buffer_lines = read_buffer_lines;    // Buffer size is dynamic
    reader->buffer = malloc(read_buffer_lines * sizeof(char *));
    reader->buffer_storage = malloc(read_buffer_lines*read_size+read_buffer_lines);
    //check mallocs
    for (int i=0; i< read_buffer_lines; i++){
        reader->buffer[i] = reader->buffer_storage + i*(read_size+1);
    }
    reader->nfiles = nfiles;              // Set the number of files}
    reader->filled = 0;                   // Initialize filled buffer slots
    // Initialize pthread mutex and condition variables
    pthread_mutex_init(&reader->mutex, NULL);
    pthread_cond_init(&reader->can_produce, NULL);
    pthread_cond_init(&reader->can_consume, NULL);
 
    reader->produce_index = 0;
    reader->consume_index = 0;
    reader->done = 0;
    return reader;
}
fastq_reader_set *  allocate_fastq_reader_set( char **barcode_filenames, char **forward_filenames, char **reverse_filenames, int nfiles, size_t read_size, size_t read_buffer_lines) {
    fastq_reader_set *reader_set=malloc(sizeof(fastq_reader_set));
    fprintf(stderr, "allocating barcode reader\n");
    reader_set->barcode_reader=allocate_fastq_reader(barcode_filenames, nfiles, 1, read_size, read_buffer_lines);
    reader_set->forward_reader=NULL;
    reader_set->reverse_reader=NULL;
    if (forward_filenames != NULL ) {
        fprintf(stderr, "allocating forward reader\n");
        reader_set->forward_reader=allocate_fastq_reader(forward_filenames, nfiles, 2, read_size, read_buffer_lines);
    }
    if (reverse_filenames != NULL) {
        fprintf(stderr, "allocating reverse reader\n");
        reader_set->reverse_reader=allocate_fastq_reader(reverse_filenames, nfiles, 3, read_size, read_buffer_lines);
    }
    return reader_set;
}
void initialize_reader_args(fastq_readers_args *reader_args, int thread_id, int set_index, int reader_type, fastq_reader_set **reader_sets) {
    reader_args->thread_id = thread_id;
    reader_args->set_index = set_index;
    reader_args->reader_type = reader_type;
    reader_args->reader_sets = reader_sets;
}
// Producer thread function
void *read_fastqs_by_set(void *arg) {
    //one set at a time
    fastq_readers_args *reader_args = (fastq_readers_args *)arg;
    const int thread_id = reader_args->thread_id;
    long long number_of_reads=0;
    fastq_reader_set **reader_sets = reader_args->reader_sets;
    const int number_of_readers = (reader_sets[0]->forward_reader != NULL) + (reader_sets[0]->reverse_reader != NULL) + 1;
    fprintf(stderr, "Thread %d reading %d files\n", thread_id, number_of_readers);
    fastq_reader *readers[number_of_readers];
    readers[0] = reader_sets[0]->barcode_reader;
    const int read_buffer_lines = readers[0]->read_buffer_lines;
    if (reader_sets[0]->forward_reader != NULL) {
        readers[1] = reader_sets[0]->forward_reader;
    }
    if (reader_sets[0]->reverse_reader != NULL) {
        readers[number_of_readers - 1] = reader_sets[0]->reverse_reader;
    }

    char line1[number_of_readers][LINE_LENGTH];
    char line2[number_of_readers][LINE_LENGTH];
    int done = 0;
    // Open the files
    for (int j=0; j<number_of_readers; j++){
        fprintf(stderr, "Opening file %s\n", readers[j]->filenames[0]);
        readers[j]->gz_pointer = gzopen(readers[j]->filenames[0], "rb");
        if (readers[j]->gz_pointer == NULL) {
            fprintf(stderr, "Error: Unable to open file %s\n", readers[j]->filenames[0]);
            exit(EXIT_FAILURE);
        }
    }
    int file_index=1;
    while(!done){
        int open_next_file=0;
        //read the lines into the given buffer
        for (int j=0; j<number_of_readers; j++){
            fastq_reader *reader=readers[j];
            if (open_next_file){
                fprintf(stderr, "1 Opening file %d %s\n", j,reader->filenames[file_index]);
                reader->gz_pointer = gzopen(reader->filenames[file_index], "rb");
                if (reader->gz_pointer == NULL) {
                    fprintf(stderr, "Error: Unable to open file %s\n", reader->filenames[file_index]);
                    exit(EXIT_FAILURE);
                } 

            }
            if ((max_reads && number_of_reads >= max_reads) || gzgets(reader->gz_pointer, line1[j], LINE_LENGTH) == Z_NULL){
                if(file_index>=reader->nfiles){
                    done=1;
                }
                else{
                    open_next_file=1;
                    number_of_reads=0;
                    fprintf(stderr, "2 Opening file %d %s\n", j,reader->filenames[file_index]);
                    reader->gz_pointer = gzopen(reader->filenames[file_index], "rb");
                    if (gzgets(reader->gz_pointer, line1[j], LINE_LENGTH) == Z_NULL){
                        perror("Error: empty FASTQ file\n");
                        exit (EXIT_FAILURE);
                    }
                }
            }
            if (!done && gzgets(reader->gz_pointer, line1[j], LINE_LENGTH) == Z_NULL){
                perror("1 Error: Incomplete record in the FASTQ file\n");
                exit (EXIT_FAILURE);
            }
            if (!done && gzgets(reader->gz_pointer, line2[j], LINE_LENGTH) == Z_NULL){
                perror("2 Error: Incomplete record in the FASTQ file\n");
                exit (EXIT_FAILURE);
            }
            if (!done && gzgets(reader->gz_pointer, line2[j], LINE_LENGTH) == Z_NULL){
                fprintf(stderr, "index %d done %d\n",j, done);
                perror("3 Error: Incomplete record in the FASTQ file\n");
                exit (EXIT_FAILURE);    
            }
            
        }
        number_of_reads++;
        if (open_next_file){
            file_index++;
        }
        // write the lines
    
        for (int j=0; j<number_of_readers; j++){
            fastq_reader *reader=readers[j];
            pthread_mutex_lock(&reader->mutex);
            if (done){
                reader->buffer[reader->produce_index][0]='\0';
                reader->buffer[(reader->produce_index+1)%read_buffer_lines][0]='\0';
                pthread_cond_signal(&reader->can_consume);
                reader->filled+=2;
                pthread_mutex_unlock(&reader->mutex);
                continue;
            }
            // Wait until there is space in the buffer
            if (reader->filled > read_buffer_lines-2) {
                pthread_cond_wait(&reader->can_produce, &reader->mutex);
            }
            // Produce lines (add it to buffer)
            strcpy(reader->buffer[reader->produce_index], line1[j]);
            strcpy(reader->buffer[(reader->produce_index+1)%read_buffer_lines], line2[j]); 
            reader->produce_index = (reader->produce_index+ 2) % read_buffer_lines;
            reader->filled+=2;
            // Signal consumer that data is available
            pthread_cond_signal(&reader->can_consume);
            pthread_mutex_unlock(&reader->mutex);
        }
    }
    
    fprintf(stderr, "Thread %d done reading\n", thread_id);
    pthread_exit(NULL);
}
void *read_fastqs(void *arg) {
    long long number_of_reads=0;
    fastq_readers_args *reader_args = (fastq_readers_args *)arg;
    const int thread_id = reader_args->thread_id;
    const int set_index = reader_args->set_index;
    const int reader_type = reader_args->reader_type;
    fastq_reader_set **reader_sets = reader_args->reader_sets;
    fastq_reader *reader=0;
    if (reader_type == 1){
        reader = reader_sets[set_index]->barcode_reader;
        DEBUG_PRINT( "Barcode reader\n");
    }
    else if (reader_type == 2){
        reader = reader_sets[set_index]->forward_reader;
        DEBUG_PRINT( "Forward reader\n");
    }
    else if (reader_type == 3){
        reader = reader_sets[set_index]->reverse_reader;
    }
    const int read_buffer_lines = reader->read_buffer_lines;

    // Open the file
    for (int i=0; i<reader->nfiles; i++){
        DEBUG_PRINT( "Opening file %s\n", reader->filenames[i]);
        reader->gz_pointer = gzopen(reader->filenames[i], "rb");
        if (reader->gz_pointer == NULL) {
            fprintf(stderr, "Error: Unable to open file %s\n", reader->filenames[i]);
            exit(EXIT_FAILURE);
        }
        char line1[LINE_LENGTH];
        char line2[LINE_LENGTH];
        int done = 0;

        
        while(!done){
            if (max_reads && number_of_reads >= max_reads){
                done=1;
            }
            if (!done && gzgets(reader->gz_pointer, line1, LINE_LENGTH) == Z_NULL){
                fprintf(stderr,"Finished reading file %s\n", reader->filenames[i]); 
                gzclose(reader->gz_pointer);
                reader->gz_pointer = NULL;
                done=1;
            }
            if (!done && gzgets(reader->gz_pointer, line1, LINE_LENGTH) == Z_NULL){
                perror("Error: Incomplete record in the FASTQ file\n");
                exit (EXIT_FAILURE);
            }
            if (!done && gzgets(reader->gz_pointer, line2, LINE_LENGTH) == Z_NULL){
                perror("Error: Incomplete record in the FASTQ file\n");
                exit (EXIT_FAILURE);
            }
            if (!done && gzgets(reader->gz_pointer, line2, LINE_LENGTH) == Z_NULL){
                perror("Error: Incomplete record in the FASTQ file\n");
                exit (EXIT_FAILURE);    
            }
            // write the lines
            pthread_mutex_lock(&reader->mutex);
            if (done){
                reader->buffer[reader->produce_index][0]='\0';
                reader->buffer[(reader->produce_index+1)%read_buffer_lines][0]='\0';
                pthread_cond_signal(&reader->can_consume);
                reader->filled+=2;
                pthread_mutex_unlock(&reader->mutex);
                break;
            }
            // Wait until there is space in the buffer
            if (reader->filled > read_buffer_lines-2) {
                pthread_cond_wait(&reader->can_produce, &reader->mutex);
            }
            // Produce lines (add it to buffer)
            strcpy(reader->buffer[reader->produce_index], line1);
            strcpy(reader->buffer[(reader->produce_index+1)%read_buffer_lines], line2); 
            reader->produce_index = (reader->produce_index+ 2) % read_buffer_lines;
            reader->filled+=2;
            number_of_reads++;
        // Signal consumer that data is available
            pthread_cond_signal(&reader->can_consume);
            pthread_mutex_unlock(&reader->mutex);
        }
    }
    reader->done = 1;
    fprintf(stderr, "Thread %d done reading\n", thread_id);
    pthread_exit(NULL);
}
void process_multiple_feature_sequences(int nsequences, char **sequences, int *orientations, feature_arrays *features, int maxHammingDistance, int nThreads, int feature_constant_offset, int max_feature_n, uint32_t *feature_index, int *hamming_distance, char *matching_sequence, uint16_t *match_position) {
    //check if the sequences are the same

    int best_hamming_distance = maxHammingDistance;
    uint32_t best_feature_index = 0;
    char best_matching_sequence[LINE_LENGTH];
    uint16_t best_match_position=0;
    char ambiguous = 0;
    int myHammingDistance=maxHammingDistance;
    for (int i=0; i<nsequences; i++){
        DEBUG_PRINT( "Sequence %s\n", sequences[i]);
        if(sequences[i]){
            uint32_t myFeatureIndex=0;
            char my_matching_sequence[LINE_LENGTH];
            uint16_t my_match_position=0;
            if (orientations[i] != 1){
                process_feature_sequence(sequences[i], features, maxHammingDistance, nThreads, feature_constant_offset, max_feature_n, &myFeatureIndex, &myHammingDistance, my_matching_sequence, &my_match_position);
                if (myFeatureIndex && myHammingDistance == best_hamming_distance) {
                    if (best_feature_index) {
                        ambiguous = 1;
                    }
                    else {
                        //first time we have a match
                        ambiguous = 0;
                        best_hamming_distance = myHammingDistance;
                        best_feature_index = myFeatureIndex;
                        best_match_position=my_match_position;
                        strcpy(best_matching_sequence, my_matching_sequence);
                    }
                }
                else if (myFeatureIndex && myHammingDistance < best_hamming_distance) {
                    ambiguous = 0;
                    best_hamming_distance = myHammingDistance;
                    best_feature_index = myFeatureIndex;
                    best_match_position=my_match_position;
                    strcpy(best_matching_sequence, my_matching_sequence);
                }
            }
            if (orientations[i] != 0){
                char reverse_sequence[LINE_LENGTH];
                strcpy(reverse_sequence, sequences[i]);
                reverse_complement_in_place(reverse_sequence);
                process_feature_sequence(reverse_sequence, features, maxHammingDistance, nThreads, feature_constant_offset, max_feature_n, &myFeatureIndex, &myHammingDistance, my_matching_sequence, &my_match_position);
            }
        }
    }
    if (ambiguous || myHammingDistance > maxHammingDistance) {
        *feature_index = 0;
        *hamming_distance = best_hamming_distance;
    } else {
        *feature_index = best_feature_index;
        *hamming_distance = best_hamming_distance;
        *match_position = best_match_position;
        if (best_feature_index) {
            strcpy(matching_sequence, best_matching_sequence);
        }
    }
}
void process_feature_sequence(char *sequence, feature_arrays *features, int maxHammingDistance, int nThreads, int feature_constant_offset, int max_feature_n, uint32_t *feature_index, int *hamming_distance, char *matching_sequence, uint16_t *match_position) {
    if (limit_search != -1 && feature_constant_offset > 0) {
        int bestHammingDistance = maxHammingDistance + 1;
        uint32_t bestFeatureIndex = 0;
        uint16_t bestMatchPosition = 0;
        char ambiguous = 0;

        for (int i = 0; i <= limit_search; i++) {
            // For i=0, this checks offset 0. For i>0, it checks +i and -i.
            for (int sign = 1; sign >= -1; sign -= 2) {
                if (i == 0 && sign == -1) continue; // Only check offset 0 once

                int offset_delta = i * sign;
                int current_offset = feature_constant_offset + offset_delta;
                if (current_offset < 0) continue;

                int currentHammingDistance = 0;
                uint32_t currentFeatureIndex = simpleCorrectFeature(sequence + current_offset, features, max_feature_n, maxHammingDistance, &currentHammingDistance);

                if (currentFeatureIndex) {
                    // If a perfect match is found, we can exit immediately.
                    if (currentHammingDistance == 0) {
                        *feature_index = currentFeatureIndex;
                        *hamming_distance = 0;
                        *match_position = current_offset;
                        memcpy(matching_sequence, sequence + current_offset, features->feature_lengths[currentFeatureIndex - 1]);
                        matching_sequence[features->feature_lengths[currentFeatureIndex - 1]] = '\0';
                        return;
                    }

                    if (currentHammingDistance < bestHammingDistance) {
                        bestHammingDistance = currentHammingDistance;
                        bestFeatureIndex = currentFeatureIndex;
                        bestMatchPosition = current_offset;
                        ambiguous = 0;
                    } else if (currentHammingDistance == bestHammingDistance) {
                        if (bestFeatureIndex != currentFeatureIndex) {
                            ambiguous = 1;
                        }
                    }
                }
            }
        }

        if (!ambiguous && bestFeatureIndex) {
            *feature_index = bestFeatureIndex;
            *hamming_distance = bestHammingDistance;
            *match_position = bestMatchPosition;
            memcpy(matching_sequence, sequence + bestMatchPosition, features->feature_lengths[bestFeatureIndex - 1]);
            matching_sequence[features->feature_lengths[bestFeatureIndex - 1]] = '\0';
        } else {
            *feature_index = 0;
            *hamming_distance = bestHammingDistance;
        }
        return;
    }
    int bestHammingDistance=0;
    int variableMaxHammingDistance=maxHammingDistance;
    uint32_t myFeatureIndex=0;
    if (feature_constant_offset > 0) {
        int constantHammingDistance=0;
        myFeatureIndex = simpleCorrectFeature(sequence + feature_constant_offset, features, max_feature_n, maxHammingDistance, &constantHammingDistance);
        if (myFeatureIndex ) {
            // need to have feature_index -1 because the feature index is 1 based but array in struct is 0 based
            memcpy(matching_sequence, sequence + feature_constant_offset, features->feature_lengths[myFeatureIndex - 1]);
            matching_sequence[features->feature_lengths[myFeatureIndex - 1]] = '\0';
            *match_position = feature_constant_offset;
            bestHammingDistance=constantHammingDistance;
            variableMaxHammingDistance=constantHammingDistance;
            //return immediately if we have a match with distance 0 (no need to check the variable hamming distance because it is the best if it also matches the offset)
            if (constantHammingDistance == 0) {
                *feature_index = myFeatureIndex;
                *hamming_distance = bestHammingDistance;
                return;
            }
        }
        //repeat with feature_constant_offset - 1 and +1 tried higher changes but seem worse
        int minusOne_feature_constant_offset = feature_constant_offset - 1;
        int minusOneFeatureIndex = simpleCorrectFeature(sequence + minusOne_feature_constant_offset, features, max_feature_n, maxHammingDistance, &constantHammingDistance);
        if (minusOneFeatureIndex && constantHammingDistance == 0) {
            *feature_index = minusOneFeatureIndex;
            *hamming_distance = bestHammingDistance;
            *match_position = minusOne_feature_constant_offset;
            return;
        }
        if (minusOneFeatureIndex && constantHammingDistance < bestHammingDistance) {
            bestHammingDistance=constantHammingDistance;
            myFeatureIndex=minusOneFeatureIndex;
            *match_position = minusOne_feature_constant_offset;
        }
        int plusOne_feature_constant_offset = feature_constant_offset + 1;
        int plusOneFeatureIndex = simpleCorrectFeature(sequence + plusOne_feature_constant_offset, features, max_feature_n, maxHammingDistance, &constantHammingDistance);
        if (plusOneFeatureIndex && constantHammingDistance == 0) {
            *feature_index = plusOneFeatureIndex;
            *hamming_distance = bestHammingDistance;
            *match_position = plusOne_feature_constant_offset;
            return;
        }
        if (plusOneFeatureIndex && constantHammingDistance < bestHammingDistance) {
            bestHammingDistance=constantHammingDistance;
            myFeatureIndex=plusOneFeatureIndex;
            *match_position = plusOne_feature_constant_offset;
        }

    }
    if (myFeatureIndex && bestHammingDistance) {
        char ambiguous = 0;
        char new_matching_sequence[LINE_LENGTH];
        uint16_t new_match_position;
        int variableHammingDistance=variableMaxHammingDistance;
        int new_feature_index = checkAndCorrectFeature(sequence, features, variableMaxHammingDistance, nThreads, &variableHammingDistance, new_matching_sequence, max_feature_n, &ambiguous, &new_match_position);
        if (variableHammingDistance < variableMaxHammingDistance) {
            myFeatureIndex = new_feature_index; // even if zero keep the new feature index because that would mean ambiguity
            strcpy(matching_sequence, new_matching_sequence);
            bestHammingDistance=variableHammingDistance;
            *match_position=new_match_position;
        }
    } else if (!myFeatureIndex) {
        char ambiguous = 0;
        char new_matching_sequence[LINE_LENGTH];
        uint16_t new_match_position;
        int variableHammingDistance=maxHammingDistance;
        int newFeatureIndex = checkAndCorrectFeature(sequence, features, variableMaxHammingDistance, nThreads, &variableHammingDistance, new_matching_sequence, max_feature_n, &ambiguous, &new_match_position);
        if (!ambiguous && newFeatureIndex && variableHammingDistance < maxHammingDistance) {
            bestHammingDistance=variableHammingDistance;
            strcpy(matching_sequence, new_matching_sequence);
            myFeatureIndex = newFeatureIndex;
            *match_position=new_match_position;
        }
    }
    if (bestHammingDistance > maxHammingDistance) {
        myFeatureIndex = 0;
    }
    *feature_index = myFeatureIndex;
    *hamming_distance = bestHammingDistance;
}
void *consume_reads(void *arg) {
    fastq_processor *processor_args = (fastq_processor *)arg;
    const int nsets = processor_args->nsets;
    fastq_reader_set **reader_sets = processor_args->reader_sets;
    const sample_args *sample_args = processor_args->sample_args;
    int done=0;
    //first one is always barcode
    char *lines_buffer=malloc(6*LINE_LENGTH);
    char *lines[6];
    char done_flags[nsets];
    memset(done_flags, 0, nsets*sizeof(char));
    const int feature_constant_offset = sample_args->feature_constant_offset;
    const int barcode_constant_offset = sample_args->barcode_constant_offset;
    feature_arrays *features = sample_args->features;
    const int maxHammingDistance = sample_args->maxHammingDistance;
    const int nThreads = sample_args->nThreads;
    int nreaders=(reader_sets[0]->forward_reader && reader_sets[0]->reverse_reader)?3:2;
    for (int i=0; i<6; i++){
        lines[i]=lines_buffer+i*LINE_LENGTH;
    }
    char **barcode_lines = lines;
    char **forward_lines = 0;
    char **reverse_lines = 0;
    if (nreaders == 3) {
        forward_lines= lines + 2;
        reverse_lines = lines + 4;
    }
    if (nreaders == 2) {
        if (reader_sets[0]->forward_reader) {
            forward_lines = lines + 2;
        }
        else {
            reverse_lines = lines + 2;
        }
    }
    while (!done) {
        // will store 0 type (barcode)
        int data_available = 0;
        for (int i = 0; i < nsets; i++) {
            if (done_flags[i]){
                continue;
            }
            fastq_reader *barcode_reader=reader_sets[i]->barcode_reader;
            fastq_reader *forward_reader=reader_sets[i]->forward_reader;
            fastq_reader *reverse_reader=reader_sets[i]->reverse_reader; 
            //check the mutexes
            if (pthread_mutex_trylock(&barcode_reader->mutex) != 0) {
                continue;
            }
            if (forward_reader && pthread_mutex_trylock(&forward_reader->mutex) != 0) {
                pthread_mutex_unlock(&barcode_reader->mutex);
                continue;
            }
            if (reverse_reader && pthread_mutex_trylock(&reverse_reader->mutex) != 0) {
                pthread_mutex_unlock(&barcode_reader->mutex);
                if (forward_reader) {
                    pthread_mutex_unlock(&forward_reader->mutex);
                }
                continue;
            }
            //check if there is data to consume
            if (barcode_reader->filled < 2) {
                pthread_mutex_unlock(&barcode_reader->mutex);
                if (forward_reader) {
                    pthread_mutex_unlock(&forward_reader->mutex);
                }
                if (reverse_reader) {
                    pthread_mutex_unlock(&reverse_reader->mutex);
                }
                continue;
            }
            if (forward_reader && forward_reader->filled < 2) {
                pthread_mutex_unlock(&barcode_reader->mutex);
                pthread_mutex_unlock(&forward_reader->mutex);
                if (reverse_reader) {
                    pthread_mutex_unlock(&reverse_reader->mutex);
                }
                continue;
            }
            if (reverse_reader && reverse_reader->filled < 2) {
                pthread_mutex_unlock(&barcode_reader->mutex);
                if (forward_reader) {
                    pthread_mutex_unlock(&forward_reader->mutex);
                }
                pthread_mutex_unlock(&reverse_reader->mutex);
                continue;
            }
            data_available = 1;
            //check if the data is null and if so set the done flag
            if (barcode_reader->buffer[barcode_reader->consume_index][0] == '\0') {
                done_flags[i] = 1;
                pthread_mutex_unlock(&barcode_reader->mutex);
                if (forward_reader) {
                    pthread_mutex_unlock(&forward_reader->mutex);
                }
                if (reverse_reader) {
                    pthread_mutex_unlock(&reverse_reader->mutex);
                }
                continue;
            }

            //copy the data to the local buffer
            strcpy(barcode_lines[0], barcode_reader->buffer[barcode_reader->consume_index]);
            strcpy(barcode_lines[1], barcode_reader->buffer[(barcode_reader->consume_index+1) % barcode_reader->read_buffer_lines]);

            if (forward_reader) {
                strcpy(forward_lines[0], forward_reader->buffer[forward_reader->consume_index]);
                strcpy(forward_lines[1], forward_reader->buffer[(forward_reader->consume_index+1) % forward_reader->read_buffer_lines]);
            }
            if (reverse_reader) {
                strcpy(reverse_lines[0], reverse_reader->buffer[reverse_reader->consume_index]);
                strcpy(reverse_lines[1], reverse_reader->buffer[(reverse_reader->consume_index+1) % reverse_reader->read_buffer_lines]);
            }
            //signal that the data has been consumed
            barcode_reader->consume_index = (barcode_reader->consume_index + 2) % barcode_reader->read_buffer_lines;
            barcode_reader->filled-=2;
            pthread_cond_signal(&barcode_reader->can_produce);
            if (forward_reader) {
                forward_reader->consume_index = (forward_reader->consume_index + 2) % forward_reader->read_buffer_lines;
                forward_reader->filled-=2;
                pthread_cond_signal(&forward_reader->can_produce);
            }
            if (reverse_reader) {
                reverse_reader->consume_index = (reverse_reader->consume_index + 2) % reverse_reader->read_buffer_lines;
                reverse_reader->filled-=2;
                pthread_cond_signal(&reverse_reader->can_produce);
            }
            pthread_mutex_unlock(&barcode_reader->mutex);
            if (forward_reader) {
                pthread_mutex_unlock(&forward_reader->mutex);
            }
            if (reverse_reader) {
                pthread_mutex_unlock(&reverse_reader->mutex);
            }
            //process the data
            char matching_sequence[LINE_LENGTH];
            int hamming_distance = 0;
            uint32_t feature_index = 0;
            uint16_t match_position = 0;
            int missing_flag=0;
            if (forward_reader && reverse_reader) {
                char *sequences[2]={0,0};
                int orientations[2]={0,0};
                sequences[0] = forward_lines[0];
                sequences[1] = reverse_lines[0];
                orientations[0] = 1;
                orientations[1] = 0;
                process_multiple_feature_sequences(2, sequences, orientations, features, maxHammingDistance, nThreads, feature_constant_offset, max_feature_n, &feature_index, &hamming_distance, matching_sequence, &match_position);
            }
            else {
                char *sequence=(forward_lines)?forward_lines[0]:reverse_lines[0];
                process_feature_sequence(sequence, features, maxHammingDistance, nThreads, feature_constant_offset, max_feature_n, &feature_index, &hamming_distance, matching_sequence, &match_position);
            }
            pthread_mutex_lock(&processor_args->process_mutex);
            if (feature_index){
                // need to have feature_index -1 because the feature index is 1 based but array in struct is 0 based
                matching_sequence[features->feature_lengths[feature_index - 1]] = '\0';               
                insert_feature_sequence(matching_sequence, feature_index, hamming_distance, match_position, sample_args->hashes, sample_args->pools);
                checkAndCorrectBarcode(barcode_lines, max_barcode_n, feature_index, match_position, sample_args->hashes, sample_args->pools, sample_args->stats, barcode_constant_offset);
            }
            else{
                missing_flag=1;
            }
            sample_args->stats->number_of_reads++;
            if (sample_args->stats->number_of_reads % 1000000 == 0)
            {
                double elapsed_time = get_time_in_seconds() - sample_args->stats->start_time;
                fprintf(stderr, "Processed %ld million reads in %.1f seconds\n", sample_args->stats->number_of_reads / 1000000, elapsed_time);
            }
            if (missing_flag){
                sample_args->stats->nMismatches++;
                sample_args->stats->total_unmatched_features++;
            }
            pthread_mutex_unlock(&processor_args->process_mutex);
        }
        for (int i=0; i<nsets; i++){
            if (!done_flags[i]){
                done=0;
                break;
            }
            done=1;
        }
        if(!done && !data_available){
            sched_yield();
        } 
    }
    //free the lines buffer
    free(lines_buffer);
    pthread_exit(NULL);
}
void free_fastq_reader(fastq_reader *reader) {
    for (int i = 0; i < reader->nfiles; i++) {
        //test if any are open
        if (reader->gz_pointer != NULL ) {
            gzclose(reader->gz_pointer);
            reader->gz_pointer = NULL;
        }
    }
    free(reader->buffer);
    free(reader->buffer_storage);
    free(reader->concatenated_filenames);
    free(reader->filenames);
    pthread_mutex_destroy(&reader->mutex);
    pthread_cond_destroy(&reader->can_produce);
    pthread_cond_destroy(&reader->can_consume);
}
void free_fastq_reader_set(fastq_reader_set *reader_set) {
    free_fastq_reader(reader_set->barcode_reader);
    free(reader_set->barcode_reader);
    if (reader_set->forward_reader) {
        free_fastq_reader(reader_set->forward_reader);
        free(reader_set->forward_reader);
    }
    if (reader_set->reverse_reader) {
        free_fastq_reader(reader_set->reverse_reader);
        free(reader_set->reverse_reader);
    }
    free(reader_set);
}
void free_fastq_processor(fastq_processor *processor_args) {
    for (int i = 0; i < processor_args->nsets; i++) {
        if (processor_args->reader_sets[i]) {
            free_fastq_reader_set(processor_args->reader_sets[i]);
        }
    }
    pthread_mutex_destroy(&processor_args->process_mutex);
}
void initiate_processor_args(fastq_processor *processor_args, sample_args *args, fastq_reader_set **reader_sets, int nsets){
    processor_args->sample_args=args;
    processor_args->reader_sets=reader_sets;
    processor_args->nsets=nsets;
    processor_args->nreaders=(reader_sets[0]->forward_reader && reader_sets[0]->reverse_reader)?3:2;
    pthread_mutex_init(&processor_args->process_mutex, NULL);
}
void process_files_in_sample(sample_args *args) {
    //allocate buffers here
    //number of lines to read into the buffer
    double  min_posterior=args->min_posterior;
    int read_by_set=(args->parallel_by_file)?0:1;
    const int sample_index = args->sample_index;
    fastq_processor processor_args;
    fastq_files_collection *fastq_files=args->fastq_files;
    const int sample_offset=fastq_files->sample_offsets[sample_index];
    const int sample_size=fastq_files->sample_sizes[sample_index];
    
    // Initialize the data structures
    fastq_reader_set* reader_sets[sample_size];
    fastq_reader_set* reader_single_set[1];
    
    fprintf(stderr, "Processing sample %d of size %d\n", sample_index, sample_size);
    if (read_by_set){
        char **forward_filenames=(fastq_files->forward_fastq)?fastq_files->forward_fastq+sample_offset:0;
        char **reverse_filenames=(fastq_files->reverse_fastq)?fastq_files->reverse_fastq+sample_offset:0;
        DEBUG_PRINT("sample size %d\n", sample_size);
        DEBUG_PRINT( "sample offset %d\n",sample_offset);
        reader_single_set[0]=allocate_fastq_reader_set(fastq_files->barcode_fastq+sample_offset,forward_filenames,reverse_filenames, sample_size, args->average_read_length, args->read_buffer_lines);
        initiate_processor_args(&processor_args, args, reader_single_set, 1);
    }
    else{
        for (int i=0; i<sample_size; i++){
            DEBUG_PRINT( "Allocating reader set %d\n", i);
            char **forward_filenames=(fastq_files->forward_fastq)?fastq_files->forward_fastq+sample_offset+i:0;
            char **reverse_filenames=(fastq_files->reverse_fastq)?fastq_files->reverse_fastq+sample_offset+i:0;
            reader_sets[i]=allocate_fastq_reader_set(fastq_files->barcode_fastq+sample_offset+i, forward_filenames, reverse_filenames, 1, args->average_read_length, args->read_buffer_lines);
        }
        initiate_processor_args(&processor_args, args, reader_sets, sample_size);
    }
    
    mkdir_p(args->directory);
        // Check if the output file directory/features_matrix.mtx exist
    const int nreaders=(fastq_files->forward_fastq && fastq_files->reverse_fastq)?3*sample_size:2*sample_size;

    

    pthread_t producer_threads[nreaders];
    int thread_index=0;

    fastq_readers_args reader_args[nreaders];
    if (read_by_set){
        initialize_reader_args(reader_args, thread_index, 0, 1, reader_single_set);
        if (pthread_create(&producer_threads[0], NULL, read_fastqs_by_set, (void *)(reader_args + thread_index)) != 0) {
            perror("Failed to create producer thread");
            exit(EXIT_FAILURE);
        }
    }
    else{
        for (int j=0; j<sample_size; j++){
            //initialize the reader args for the barcode reader
            fprintf(stderr, "Creating barcode reader thread %d set index %d\n", thread_index, j);
            initialize_reader_args(reader_args+thread_index, thread_index, j, 1, reader_sets);
            if (pthread_create(&producer_threads[thread_index], NULL, read_fastqs, (void *)(reader_args + thread_index)) != 0) {
                perror("Failed to create producer thread");
                exit(EXIT_FAILURE);
            }
            thread_index++;
            if (fastq_files->forward_fastq){
                fprintf(stderr, "Creating forward reader thread %d set index %d\n", thread_index, j);
                initialize_reader_args(reader_args+thread_index, thread_index, j, 2, reader_sets);
                if (pthread_create(&producer_threads[thread_index], NULL, read_fastqs, (void *)(reader_args + thread_index)) != 0) {
                    perror("Failed to create producer thread");
                    exit(EXIT_FAILURE);
                }
                thread_index++;
            }
            if (fastq_files->reverse_fastq){
                initialize_reader_args(reader_args+thread_index, thread_index, j, 3, reader_sets);
                if (pthread_create(&producer_threads[thread_index], NULL, read_fastqs, (void *)(reader_args + thread_index)) != 0) {
                    perror("Failed to create producer thread");
                    exit(EXIT_FAILURE);
                }
                thread_index++;
            }      
        }
    }
    const int nconsumers=args->consumer_threads_per_set;
    pthread_t consumer_threads[nconsumers];
    fprintf(stderr, "Will use %d threads to process the reads\n", nconsumers);
    for (int j=0; j<nconsumers; j++){
        if (pthread_create(&consumer_threads[j], NULL, consume_reads, (void *)&processor_args) != 0) {
            perror("Failed to create consumer thread");
            exit(EXIT_FAILURE);
        }
    }
    //join the threads
    const int nproducer_threads=(read_by_set)?1:nreaders;
    for (int j=0; j<nproducer_threads; j++){
        pthread_join(producer_threads[j], NULL);
    }
    for (int j=0; j<nconsumers; j++){
        pthread_join(consumer_threads[j], NULL);
    }
    fprintf(stderr, "Finished processing sample %d\n", sample_index);
    finalize_processing(args->features, args->hashes, args->directory, args->pools,args->stats, args->stringency, args->min_counts, min_posterior);
    free_fastq_processor(&processor_args);
}
void initialize_data_structures(data_structures *hashes){
    hashes->filtered_hash = g_hash_table_new(hash_int32, equal_int32);
    hashes->unique_features_match = g_hash_table_new(g_str_hash, g_str_equal);    
    // Use g_hash_table_new_full to provide our custom value destroyer function.
    hashes->sequence_umi_hash = g_hash_table_new_full(hash_int64, equal_int64, NULL, destroy_feature_umi_counts);
    
    hashes->neighbors_queue=malloc(sizeof(Queue));
    //check if the memory allocation was successful
    if (hashes->neighbors_queue == NULL) {
        perror("Failed to allocate memory for neighbors queue");
        exit(EXIT_FAILURE);
    }
    init_queue(hashes->neighbors_queue);
}
void destroy_data_structures(data_structures *hashes){
    g_hash_table_destroy(hashes->filtered_hash);
    g_hash_table_destroy(hashes->unique_features_match);
    g_hash_table_destroy(hashes->sequence_umi_hash);
    free_queue(hashes->neighbors_queue);
    free(hashes->neighbors_queue);
}


int existing_output_skip(char keep_existing, char *directory){
    if (keep_existing && file_exists(directory)){
        char matrix_filename[4096];
        strcpy(matrix_filename, directory);
        strcat(matrix_filename, "features_matrix.mtx");
        if (file_exists(matrix_filename)){
            fprintf(stderr, "Matrix file %s found and skipping %s\n", matrix_filename, get_basename(directory));
            return 1;
        }
    }
    return 0;
}
size_t get_file_size(char *filepath){
    char *filename=filepath;
    struct stat file_stat;
    if (stat(filename, &file_stat) == 0) {
        return file_stat.st_size;
    } else {
        perror("stat");
    }
    return 0;
}

int compare_file_sizes(const void *a, const void *b, void *context) {
    size_t *sizes = (size_t *)context;
    int index_a = *(int *)a;
    int index_b = *(int *)b;

    // Compare sizes in descending order (so larger sizes come first)
    if (sizes[index_a] < sizes[index_b]) {
        return 1;
    } else if (sizes[index_a] > sizes[index_b]) {
        return -1;
    } else {
        return 0;
    }
}

// Custom qsort that allows passing a context (like the sizes array)
void qsort_with_context(void *base, size_t num, size_t size, 
                        int (*compar)(const void *, const void *, void *), 
                        void *context) {
    char *ptr = (char *)base;
    for (size_t i = 0; i < num - 1; i++) {
        for (size_t j = i + 1; j < num; j++) {
            if (compar(ptr + i * size, ptr + j * size, context) > 0) {
                // Swap the elements
                char temp[size];
                memcpy(temp, ptr + i * size, size);
                memcpy(ptr + i * size, ptr + j * size, size);
                memcpy(ptr + j * size, temp, size);
            }
        }
    }
}

// Main function to sort the indices of samples by their size and store the order in sample_order
void sort_samples_by_size(fastq_files_collection *fastq_files, int *sample_order) {
    // Allocate memory for the sizes array to store the total size of each sample
    size_t *sizes = malloc(fastq_files->nsamples * sizeof(size_t));
    if (sizes == NULL) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }

    // Calculate the total size for each sample and store it in the sizes array
    int sample_index = 0;
    for (int i = 0; i < fastq_files->nsamples; i++) {
        sizes[i] = 0;  // Initialize size for each sample
        for (int j = 0; j < fastq_files->sample_sizes[i]; j++) {
            sizes[i] += get_file_size(fastq_files->barcode_fastq[sample_index]);
            if (fastq_files->forward_fastq) {
                sizes[i] += get_file_size(fastq_files->forward_fastq[sample_index]);
            }
            if (fastq_files->reverse_fastq) {
                sizes[i] += get_file_size(fastq_files->reverse_fastq[sample_index]);
            }
            sample_index++;
        }
        sample_order[i] = i;  // Initialize sample_order with indices 0 to nsamples-1
    }

    // Sort the sample_order array based on the sizes array
    qsort_with_context(sample_order, fastq_files->nsamples, sizeof(int), compare_file_sizes, sizes);

    // Free the allocated memory for sizes
    free(sizes);
}
void cleanup_sample(memory_pool_collection *pools, data_structures *hashes){
    destroy_data_structures(hashes);      /* first: walk and free hash nodes */
    free_memory_pool_collection(pools);   /* then: free the underlying pools */
}
int calculate_threads(fastq_files_collection *fastq_files, int parallel_by_file, int index, int threads_used, int available_threads, int *consumer_threads_per_set, int *search_threads_per_consumer, int *max_concurrent_processes, int set_consumer_threads_per_set, int set_search_threads_per_consumer){
    //assume that all concurrent processes will eventually become available for calculations
    int nsamples=fastq_files->nsamples-index;
    int sample_index=fastq_files->sorted_index[index];
    int sample_size=fastq_files->sample_sizes[sample_index];
    int producer_threads=(parallel_by_file)?sample_size:1;
    int remaining_threads=available_threads-producer_threads-threads_used;
    fprintf(stderr, "Remaining threads %d\n", remaining_threads);
    double best_value=0;
    int best_search=4; //default to 1 in worst case
    int best_consumer=1;
    int best_concurrent=1;
    double search_value[3]={3.2,1.9,1.0};
    int sthreads[3]={4,2,1};
    double producer_cost=0.6*(double)producer_threads;
    int start_sthread_index=0;
    int end_sthread_index=2;
    int start_consumer=1;
    int end_consumer=(remaining_threads)?remaining_threads:1;
    if (set_consumer_threads_per_set){
        start_consumer=set_consumer_threads_per_set;
        end_consumer=set_consumer_threads_per_set;
    }
    //set the search threads if specified by user
    if (set_search_threads_per_consumer){
        if (set_search_threads_per_consumer<=1){
            start_sthread_index=2;
            end_sthread_index=2;
        }
        else if (set_search_threads_per_consumer<=3){
            start_sthread_index=1;
            end_sthread_index=1;
        }
        else {
            start_sthread_index=0;
            end_sthread_index=0;
        }
    }    
    for (int i=start_sthread_index; i<=end_sthread_index; i++){   
        for (int consumer=start_consumer; consumer<=end_consumer; consumer++){
            for (int concurrent=1; concurrent<=*max_concurrent_processes && concurrent <= nsamples; concurrent++){
                double dconsumer=(double) consumer;
                double dsearch=(double)search_value[i];
                double dconcurrent=(double)concurrent;
                double value=dconsumer*dsearch*dconcurrent-dconsumer*dconsumer/4.0;
               
                DEBUG_PRINT( "consumer %d search %d concurrent %d Value %f\n",consumer,sthreads[i],concurrent,value);
              
                if (value>best_value){
                    best_value=value;
                    best_search=sthreads[i];
                    best_consumer=consumer;
                    best_concurrent=concurrent;
                    DEBUG_PRINT( "Best value %f consumer %d search %d concurrent %d\n", best_value, best_consumer, best_search, best_concurrent);
                }
            }
        }
    }
    *consumer_threads_per_set=best_consumer;
    *search_threads_per_consumer=best_search;
    *max_concurrent_processes=best_concurrent;
    DEBUG_PRINT( "Recalculating: Use %d consumer threads, %d search threads\n", best_consumer, best_search );
    return best_consumer*best_search + producer_cost;
}

int calculate_initial_threads(fastq_files_collection *fastq_files, int parallel_by_file, int available_threads, int *consumer_threads_per_set, int *search_threads_per_consumer, int *max_concurrent_processes, int set_consumer_threads_per_set, int set_search_threads_per_consumer){
    int nsamples=fastq_files->nsamples;
    int sample_size=fastq_files->max_sample_size;
    //start off with producer threads
    int producer_threads=(parallel_by_file)?sample_size:1;
    int remaining_threads=available_threads-producer_threads; 
    //search threads should be 4:2:1 depending on the number of consumer threads available - could use 3 but that is a waste
    //metric for work done is approximately consumer*search
    double best_value=0;
    int best_search=4; //default to 1 in worst case
    int best_consumer=1;
    int best_concurrent=1;
    //int search=(remaining_threads>=4)?4:2;
    //if (remaining_threads<2){
    //    search=1;
    //}
    //prefer high concurrent and low consumer - less overhead - put consumer in outer loop

    double search_value[3]={3.2,1.9,1.0};
    int sthreads[3]={4,2,1};
    double producer_cost=0.6*(double)producer_threads;
    int start_sthread_index=0;
    int end_sthread_index=2;
    int start_consumer=1;
    int end_consumer=(remaining_threads)?remaining_threads:1;
    if (set_consumer_threads_per_set){
        start_consumer=set_consumer_threads_per_set;
        end_consumer=set_consumer_threads_per_set;
    }
    //set the search threads if specified by user
    if (set_search_threads_per_consumer){
        if (set_search_threads_per_consumer<=1){
            start_sthread_index=2;
            end_sthread_index=2;
        }
        else if (set_search_threads_per_consumer<=3){
            start_sthread_index=1;
            end_sthread_index=1;
        }
        else {
            start_sthread_index=0;
            end_sthread_index=0;
        }
    }    
    for (int i=start_sthread_index; i<=end_sthread_index; i++){   
        for (int consumer=start_consumer; consumer<=end_consumer; consumer++){
            for (int concurrent=1; concurrent<=*max_concurrent_processes && concurrent <= nsamples; concurrent++){
                //check that this is less than number of threads - remaining threads instead of available for same reason
                if (((double)consumer*sthreads[i]+producer_cost)*(double)concurrent>available_threads){
                    break;
                }
                double dconsumer=(double) consumer;
                double dsearch=(double)search_value[i];
                double dconcurrent=(double)concurrent;
                double value=dconsumer*dsearch*dconcurrent-dconsumer*dconsumer/4.0;
                DEBUG_PRINT( "consumer %d search %d concurrent %d Value %f\n",consumer,sthreads[i],concurrent,value);

                if (value>best_value){
                    best_value=value;
                    best_search=sthreads[i];
                    best_consumer=consumer;
                    best_concurrent=concurrent;
                    DEBUG_PRINT( "Best value %f consumer %d search %d concurrent %d\n", best_value, best_consumer, best_search, best_concurrent);

                }
            }
        }
    }
    *consumer_threads_per_set=best_consumer;
    *search_threads_per_consumer=best_search;
    *max_concurrent_processes=best_concurrent;
    return best_consumer*best_search + producer_cost;
}
void reverse_in_place(char *str) {
    int left = 0;
    int right = strlen(str) - 1;
    // Swap characters from both ends of the string
    while (left < right) {
        char temp = str[left];
        str[left] = str[right];
        str[right] = temp;

        left++;
        right--;
    }
}
char complement(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return base;  // For unexpected characters, return the base itself
    }
}

// Function to calculate the reverse complement of a sequence in place
void reverse_complement_in_place(char *seq) {
    int left = 0;
    int right = strlen(seq) - 1;

    // Loop through the sequence from both ends, swapping and complementing
    while (left <= right) {
        // Get the complement of both ends
        char left_complement = complement(seq[left]);
        char right_complement = complement(seq[right]);

        // Swap the left and right complements
        seq[left] = right_complement;
        seq[right] = left_complement;

        left++;
        right--;
    }
}

void reverse_complement_sequence(char *sequence,  char *reverse, int length){
    for (int i = 0; i < length; i++) {
        const char comp_base =match[(unsigned char)sequence[length - i - 1]];
        if (comp_base == 0){
            fprintf(stderr, "Error: Invalid base %c in sequence %s\n", sequence[length - i - 1], sequence);
            exit(EXIT_FAILURE);
        }
        reverse[i] = comp_base;
    }
    reverse[length] = '\0';
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
void populate_sample_args(sample_args *args, int sample_index,char *directory, fastq_files_collection *fastq_files, feature_arrays *features, int maxHammingDistance, int nThreads, memory_pool_collection *pools, statistics *stats, data_structures *hashes, uint16_t stringency, uint16_t min_counts, int barcode_constant_offset, int feature_constant_offset, int read_buffer_lines, int average_read_length, int parallel_by_file, double min_posterior, int consumer_threads_per_set){ 
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
                args->parallel_by_file=parallel_by_file;
                args->min_posterior=min_posterior;
                args->consumer_threads_per_set=consumer_threads_per_set;

            }

