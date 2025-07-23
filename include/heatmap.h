#ifndef HEATMAP_H
#define HEATMAP_H

#include "common.h"

void generate_heatmap(const char *directory, feature_arrays *features, int **coexpression_histograms);
void generate_deduped_heatmap(const char *directory, feature_arrays *features, GArray **feature_hist, int *total_deduped_counts, int histogram_minimum_counts);

#endif /* HEATMAP_H */