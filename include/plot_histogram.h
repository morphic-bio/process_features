#ifndef PLOT_HISTOGRAM_H
#define PLOT_HISTOGRAM_H

#include "common.h"
#include <glib.h>
#include "EMfit.h" // For NBSignalCut

// Function to plot combined cumulative/average dedupe histogram with EM fit and cutoffs
void plot_combined_histogram_with_em(const char *directory,
                                     GArray **feature_hist,
                                     NBSignalCut cumulative_em_fit,
                                     double posterior_cutoff,
                                     int n_features,
                                     double em_cumulative_limit);

// Helper function to calculate theoretical fit values
void calculate_theoretical_fit(NBSignalCut em_fit, 
                              int max_count, 
                              double *fit_values);

void plot_simple_histogram(const char *directory,
                           const char *filename_leaf,
                           const char *title,
                           const char *xtitle,
                           const char *ytitle,
                           GArray *hist_data);

#endif // PLOT_HISTOGRAM_H
