#ifndef HEATMAP_H
#define HEATMAP_H

#include "common.h"

void generate_heatmap(const char *directory, feature_arrays *features, int **coexpression_histograms);

#endif /* HEATMAP_H */