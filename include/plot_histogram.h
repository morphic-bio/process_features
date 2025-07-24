#ifndef PLOT_HISTOGRAM_H
#define PLOT_HISTOGRAM_H

#include "common.h"

// Function to plot combined cumulative/average dedupe histogram with EM fit and cutoffs
void plot_combined_histogram_with_em(const char *directory,
                                     GArray *histogram,
                                     NBSignalCut em_fit,
                                     uint16_t min_counts,
                                     double posterior_cutoff,
                                     int n_features,
                                     double em_cumulative_limit);

// Helper function to calculate theoretical fit values
void calculate_theoretical_fit(NBSignalCut em_fit, 
                              int max_count, 
                              double *fit_values);

#endif // PLOT_HISTOGRAM_H
