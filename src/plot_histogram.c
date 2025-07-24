#include "../include/plot_histogram.h"
#include "../include/globals.h"
#include "../include/EMfit.h"
#include <math.h>

// Helper function to calculate negative binomial PMF
double nb_pmf(int k, double r, double p) {
    if (k < 0 || r <= 0.0 || p <= 0.0 || p >= 1.0) return 0.0;
    return exp(lgamma(k + r) - lgamma(k + 1) - lgamma(r) + r*log(p) + k*log(1.0 - p));
}

// Helper function to calculate Poisson PMF
double poisson_pmf(int k, double lambda) {
    if (k < 0 || lambda <= 0) return 0.0;
    return exp(k*log(lambda) - lambda - lgamma(k + 1));
}

void calculate_theoretical_fit(NBSignalCut em_fit, int max_count, double *fit_values) {
    for (int k = 0; k <= max_count; k++) {
        fit_values[k] = 0.0;
        
        for (int comp = 0; comp < em_fit.n_comp; comp++) {
            double component_value = 0.0;
            
            if (em_fit.p[comp] == -1.0) {
                // Poisson component (sentinel value -1.0 for p)
                component_value = poisson_pmf(k, em_fit.r[comp]);
            } else {
                // Negative binomial component
                component_value = nb_pmf(k, em_fit.r[comp], em_fit.p[comp]);
            }
            
            fit_values[k] += em_fit.weight[comp] * component_value;
        }
    }
}

// New function to calculate individual component values
void calculate_individual_components(NBSignalCut em_fit, int max_count, double **component_values, long total_counts) {
    for (int comp = 0; comp < em_fit.n_comp; comp++) {
        for (int k = 0; k <= max_count; k++) {
            double component_value = 0.0;
            
            if (em_fit.p[comp] == -1.0) {
                // Poisson component
                component_value = poisson_pmf(k, em_fit.r[comp]);
            } else {
                // Negative binomial component
                component_value = nb_pmf(k, em_fit.r[comp], em_fit.p[comp]);
            }
            
            // Scale by weight and total counts
            component_values[comp][k] = em_fit.weight[comp] * component_value * total_counts;
        }
    }
}

void generate_plotly_html(const char *filename, 
                         const uint32_t *hist_data, 
                         int hist_len,
                         NBSignalCut em_fit, 
                         uint16_t min_counts,
                         double posterior_cutoff) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot create file %s\n", filename);
        return;
    }

    // Calculate theoretical fit values
    double *fit_values = malloc((hist_len + 1) * sizeof(double));
    if (!fit_values) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        fclose(fp);
        return;
    }
    
    // Allocate memory for individual components
    double **component_values = malloc(em_fit.n_comp * sizeof(double*));
    for (int comp = 0; comp < em_fit.n_comp; comp++) {
        component_values[comp] = malloc((hist_len + 1) * sizeof(double));
    }
    
    calculate_theoretical_fit(em_fit, hist_len, fit_values);
    
    // Find total counts for scaling the theoretical fit
    long total_counts = 0;
    for (int i = 0; i < hist_len; i++) {
        total_counts += hist_data[i];
    }

    // Calculate individual component values
    calculate_individual_components(em_fit, hist_len, component_values, total_counts);

    // Write HTML header
    fprintf(fp, "<!DOCTYPE html>\n");
    fprintf(fp, "<html>\n");
    fprintf(fp, "<head>\n");
    fprintf(fp, "    <title>Average Feature UMI counts</title>\n");
    fprintf(fp, "    <script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script>\n");
    fprintf(fp, "</head>\n");
    fprintf(fp, "<body>\n");
    fprintf(fp, "    <div id=\"plot\" style=\"width:100%%;height:600px;\"></div>\n");
    fprintf(fp, "    <script>\n");

    // Histogram data
    fprintf(fp, "        var hist_x = [");
    for (int i = 0; i < hist_len; i++) {
        if (i > 0) fprintf(fp, ", ");
        fprintf(fp, "%d", i);
    }
    fprintf(fp, "];\n");

    fprintf(fp, "        var hist_y = [");
    for (int i = 0; i < hist_len; i++) {
        if (i > 0) fprintf(fp, ", ");
        fprintf(fp, "%u", hist_data[i]);
    }
    fprintf(fp, "];\n");

    // Calculate y-axis range for plots
    fprintf(fp, "        var y_axis_max = Math.max(...hist_y) * 1.1;\n");
    fprintf(fp, "        var log_y_min = Math.log10(0.5);\n");
    fprintf(fp, "        var log_y_max = (y_axis_max > 0) ? Math.log10(y_axis_max) : 0;\n");


    // Find last bin with frequency >= 1 to set x-axis range
    int last_bin_with_data = 0;
    for (int i = hist_len - 1; i >= 1; i--) {
        if (hist_data[i] >= 1) {
            last_bin_with_data = i;
            break;
        }
    }
    int x_axis_max = (last_bin_with_data + 4) / 5 * 5;
    if (x_axis_max == 0) x_axis_max = 5;


    // Theoretical fit data (scaled to match histogram counts)
    fprintf(fp, "        var fit_x = [");
    for (int i = 0; i < hist_len; i++) {
        if (i > 0) fprintf(fp, ", ");
        fprintf(fp, "%d", i);
    }
    fprintf(fp, "];\n");

    fprintf(fp, "        var fit_y = [");
    for (int i = 0; i < hist_len; i++) {
        if (i > 0) fprintf(fp, ", ");
        fprintf(fp, "%.6f", fit_values[i] * total_counts);
    }
    fprintf(fp, "];\n");

    // Individual component data
    const char* component_colors[] = {"blue", "green", "orange", "purple"};
    const char* component_dash[] = {"dot", "dashdot", "longdash", "longdashdot"};
    
    for (int comp = 0; comp < em_fit.n_comp; comp++) {
        fprintf(fp, "        var comp%d_x = [", comp);
        for (int i = 0; i < hist_len; i++) {
            if (i > 0) fprintf(fp, ", ");
            fprintf(fp, "%d", i);
        }
        fprintf(fp, "];\n");

        fprintf(fp, "        var comp%d_y = [", comp);
        for (int i = 0; i < hist_len; i++) {
            if (i > 0) fprintf(fp, ", ");
            fprintf(fp, "%.6f", component_values[comp][i]);
        }
        fprintf(fp, "];\n");
    }

    // Define traces
    fprintf(fp, "        var histogram_trace = {\n");
    fprintf(fp, "            x: hist_x,\n");
    fprintf(fp, "            y: hist_y,\n");
    fprintf(fp, "            type: 'bar',\n");
    fprintf(fp, "            name: 'Observed Histogram',\n");
    fprintf(fp, "            marker: { color: 'rgba(55, 128, 191, 0.7)', line: { color: 'rgba(55, 128, 191, 1.0)', width: 1 } },\n");
    fprintf(fp, "            opacity: 0.7\n");
    fprintf(fp, "        };\n");

    fprintf(fp, "        var fit_trace = {\n");
    fprintf(fp, "            x: fit_x,\n");
    fprintf(fp, "            y: fit_y,\n");
    fprintf(fp, "            type: 'scatter',\n");
    fprintf(fp, "            mode: 'lines',\n");
    fprintf(fp, "            name: 'Overall EM Fit (%d components)',\n", em_fit.n_comp);
    fprintf(fp, "            line: { color: 'red', width: 3 }\n");
    fprintf(fp, "        };\n");

    // Individual component traces
    for (int comp = 0; comp < em_fit.n_comp; comp++) {
        //const char* dist_type = (em_fit.p[comp] == -1.0) ? "Poisson" : "NB";
        const char* component_name = "";
        
        // Determine component name based on distribution and parameters
        if (em_fit.p[comp] == -1.0) {
            component_name = "Background (Poisson)";
        } else {
            // For NB components, try to identify based on mean
            double mean = em_fit.r[comp] * (1.0 - em_fit.p[comp]) / em_fit.p[comp];
            if (comp == 0 || mean < 10) {
                component_name = "Signal (NB)";
            } else {
                component_name = "Multiplet (NB)";
            }
        }
        
        fprintf(fp, "        var comp%d_trace = {\n", comp);
        fprintf(fp, "            x: comp%d_x,\n", comp);
        fprintf(fp, "            y: comp%d_y,\n", comp);
        fprintf(fp, "            type: 'scatter',\n");
        fprintf(fp, "            mode: 'lines',\n");
        if (em_fit.p[comp] == -1.0) {
            fprintf(fp, "            name: 'Component %d: %s (λ=%.3f, w=%.3f)',\n", comp + 1, component_name, em_fit.r[comp], em_fit.weight[comp]);
        } else {
            fprintf(fp, "            name: 'Component %d: %s (r=%.3f, p=%.3f, w=%.3f)',\n", comp + 1, component_name, em_fit.r[comp], em_fit.p[comp], em_fit.weight[comp]);
        }
        fprintf(fp, "            line: { color: '%s', width: 2, dash: '%s' },\n", component_colors[comp % 4], component_dash[comp % 4]);
        fprintf(fp, "            opacity: 0.8\n");
        fprintf(fp, "        };\n");
    }

    // Vertical lines for cutoffs
    fprintf(fp, "        var min_cutoff_trace = {\n");
    fprintf(fp, "            x: [%d, %d],\n", em_fit.k_min_signal, em_fit.k_min_signal);
    fprintf(fp, "            y: [0, y_axis_max],\n");
    fprintf(fp, "            type: 'scatter',\n");
    fprintf(fp, "            mode: 'lines',\n");
    fprintf(fp, "            name: 'Min Signal Cutoff (%d)',\n", em_fit.k_min_signal);
    fprintf(fp, "            line: { color: 'darkgreen', width: 2, dash: 'dash' },\n");
    fprintf(fp, "            showlegend: true\n");
    fprintf(fp, "        };\n");

    fprintf(fp, "        var max_cutoff_trace = {\n");
    fprintf(fp, "            x: [%d, %d],\n", em_fit.k_max_signal, em_fit.k_max_signal);
    fprintf(fp, "            y: [0, y_axis_max],\n");
    fprintf(fp, "            type: 'scatter',\n");
    fprintf(fp, "            mode: 'lines',\n");
    fprintf(fp, "            name: 'Max Signal Cutoff (%d)',\n", em_fit.k_max_signal);
    fprintf(fp, "            line: { color: 'darkorange', width: 2, dash: 'dash' },\n");
    fprintf(fp, "            showlegend: true\n");
    fprintf(fp, "        };\n");

    // Dummy trace for posterior cutoff
    fprintf(fp, "        var posterior_trace = {\n");
    fprintf(fp, "            x: [null],\n");
    fprintf(fp, "            y: [null],\n");
    fprintf(fp, "            type: 'scatter',\n");
    fprintf(fp, "            mode: 'markers',\n");
    fprintf(fp, "            name: 'Posterior Cutoff: %.3f',\n", posterior_cutoff);
    fprintf(fp, "            marker: { size: 0 },\n");
    fprintf(fp, "            showlegend: true\n");
    fprintf(fp, "        };\n");

    // Layout
    fprintf(fp, "        var layout = {\n");
    fprintf(fp, "            title: {\n");
    fprintf(fp, "                text: 'Average Feature UMI counts<br><sub>BIC: %.2f, Components: %d</sub>',\n", em_fit.bic, em_fit.n_comp);
    fprintf(fp, "                x: 0.5\n");
    fprintf(fp, "            },\n");
    fprintf(fp, "            xaxis: {\n");
    fprintf(fp, "                title: 'UMI counts',\n");
    fprintf(fp, "                type: 'linear',\n");
    fprintf(fp, "                range: [0.5, %d]\n", x_axis_max);
    fprintf(fp, "            },\n");
    fprintf(fp, "            yaxis: {\n");
    fprintf(fp, "                title: 'Frequency',\n");
    fprintf(fp, "                type: 'linear',\n");
    fprintf(fp, "                autorange: true\n");
    fprintf(fp, "            },\n");
    fprintf(fp, "            updatemenus: [\n");
    fprintf(fp, "                {\n");
    fprintf(fp, "                    buttons: [\n");
    fprintf(fp, "                        {\n");
    fprintf(fp, "                            args: [{'yaxis.type': 'linear', 'yaxis.autorange': true}],\n");
    fprintf(fp, "                            label: 'Linear Y-axis',\n");
    fprintf(fp, "                            method: 'relayout'\n");
    fprintf(fp, "                        },\n");
    fprintf(fp, "                        {\n");
    fprintf(fp, "                            args: [{'yaxis.type': 'log', 'yaxis.range': [log_y_min, log_y_max]}],\n");
    fprintf(fp, "                            label: 'Log Y-axis',\n");
    fprintf(fp, "                            method: 'relayout'\n");
    fprintf(fp, "                        }\n");
    fprintf(fp, "                    ],\n");
    fprintf(fp, "                    direction: 'down',\n");
    fprintf(fp, "                    pad: {t: 10, r: 10},\n");
    fprintf(fp, "                    showactive: true,\n");
    fprintf(fp, "                    type: 'dropdown',\n");
    fprintf(fp, "                    x: 0.1,\n");
    fprintf(fp, "                    xanchor: 'left',\n");
    fprintf(fp, "                    y: 1.15,\n");
    fprintf(fp, "                    yanchor: 'top'\n");
    fprintf(fp, "                }\n");
    fprintf(fp, "            ],\n");
    fprintf(fp, "            showlegend: true,\n");
    fprintf(fp, "            legend: {\n");
    fprintf(fp, "                x: 1.02,\n");
    fprintf(fp, "                y: 1.0,\n");
    fprintf(fp, "                xanchor: 'left',\n");
    fprintf(fp, "                yanchor: 'top',\n");
    fprintf(fp, "                bgcolor: 'rgba(255,255,255,0.8)',\n");
    fprintf(fp, "                bordercolor: 'black',\n");
    fprintf(fp, "                borderwidth: 1\n");
    fprintf(fp, "            },\n");
    fprintf(fp, "            hovermode: 'x unified'\n");
    fprintf(fp, "        };\n");

    // Combine all traces
    fprintf(fp, "        var data = [histogram_trace, fit_trace");
    for (int comp = 0; comp < em_fit.n_comp; comp++) {
        fprintf(fp, ", comp%d_trace", comp);
    }
    fprintf(fp, ", min_cutoff_trace, max_cutoff_trace, posterior_trace];\n");
    
    fprintf(fp, "        Plotly.newPlot('plot', data, layout);\n");

    fprintf(fp, "    </script>\n");
    fprintf(fp, "</body>\n");
    fprintf(fp, "</html>\n");

    fclose(fp);
    
    // Free allocated memory
    free(fit_values);
    for (int comp = 0; comp < em_fit.n_comp; comp++) {
        free(component_values[comp]);
    }
    free(component_values);
    
    fprintf(stderr, "Average histogram plot with individual components saved to: %s\n", filename);
}

void plot_average_histogram_with_em(const char *directory,
                                      GArray *histogram,
                                      NBSignalCut em_fit,
                                      uint16_t min_counts,
                                      double posterior_cutoff,
                                      int n_features,
                                      double em_cumulative_limit) {
    if (!histogram || histogram->len == 0) {
        fprintf(stderr, "Warning: Empty histogram provided for plotting\n");
        return;
    }
    if (n_features == 0) {
        fprintf(stderr, "Warning: Number of features is 0, cannot generate average histogram.\n");
        return;
    }

    // Create a new GArray for the average histogram
    uint32_t *avg_hist_data = malloc(histogram->len * sizeof(uint32_t));
    if (!avg_hist_data) {
        fprintf(stderr, "Error: Memory allocation for average histogram failed.\n");
        return;
    }
    
    long total_avg_counts = 0;
    for (guint i = 0; i < histogram->len; i++) {
        avg_hist_data[i] = (uint32_t)round((double)g_array_index(histogram, uint32_t, i) / n_features);
        total_avg_counts += avg_hist_data[i];
    }
    
    // Recalculate cutoffs using the average histogram
    NBSignalCut avg_fit = em_fit; // Copy the model parameters
    avg_fit.total_counts_in_hist = total_avg_counts;
    determine_signal_cutoff_from_fit(&avg_fit, histogram->len, posterior_cutoff, min_counts, em_cumulative_limit);

    // Create output filename
    char filename[FILENAME_LENGTH];
    sprintf(filename, "%s/average_histogram_with_em_fit.html", directory);

    // Generate the plot using the average data and recalculated cutoffs
    generate_plotly_html(filename, avg_hist_data, histogram->len, avg_fit, min_counts, posterior_cutoff);
    
    // Free the allocated memory for the average histogram
    free(avg_hist_data);
}
