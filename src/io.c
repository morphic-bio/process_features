#include "../include/common.h"
#include "../include/globals.h"
#include "../include/prototypes.h"
#include "../include/utils.h"
#include "../include/io.h"

// All non-I/O functions from the original assignBarcodes.c remain here.
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

void free_memory_pool_storage(memory_pool *pool) {
    storage_block *current = pool->first_block;
    while (current != NULL) {
        storage_block *next = current->next;
        free(current->storage);
        free(current);
        current = next;
    }
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

// ... other non-I/O functions