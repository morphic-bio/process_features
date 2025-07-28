#include "../include/plot_histogram.h"
#include "../include/globals.h"
#include "../include/EMfit.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>


// --- Static Helper Functions for Plotting ---

// Helper function to calculate negative binomial PMF
static double nb_pmf(int k, double r, double p) {
    if (k < 0 || r <= 0.0 || p <= 0.0 || p >= 1.0) return 0.0;
    return exp(lgamma(k + r) - lgamma(k + 1) - lgamma(r) + r*log(p) + k*log(1.0 - p));
}

// Helper function to calculate Poisson PMF
static double poisson_pmf(int k, double lambda) {
    if (k < 0 || lambda <= 0) return 0.0;
    return exp(k*log(lambda) - lambda - lgamma(k + 1));
}

// Function to calculate individual component values
static void calculate_individual_components(NBSignalCut em_fit, int max_count, double **component_values, long total_counts) {
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

static void print_js_array_uint32(FILE* f, const char* name, const uint32_t* data, int len) {
    fprintf(f, "var %s = [", name);
    for(int i=0; i<len; ++i) {
        fprintf(f, "%u%s", data[i], (i==len-1 ? "" : ","));
    }
    fprintf(f, "];\n");
}

static void print_js_array_double(FILE* f, const char* name, const double* data, int len) {
    fprintf(f, "var %s = [", name);
    for(int i=0; i<len; ++i) {
        fprintf(f, "%.6f%s", data[i], (i==len-1 ? "" : ","));
    }
    fprintf(f, "];\n");
}

static void generate_traces_for_plotly(FILE *fp, const char* prefix, const uint32_t* hist_data, int hist_len, NBSignalCut fit) {
    long total_counts = 0;
    for (int i = 0; i < hist_len; i++) total_counts += hist_data[i];
    
    const char* view_name = (strcmp(prefix, "cumulative") == 0) ? "Cumulative" : "Average";

    // --- Histogram Trace ---
    fprintf(fp, "var %s_hist_trace = { x: x_axis, y: %s_hist_y, type: 'bar', name: 'Observed Histogram', legendgroup: '%s', marker: { color: 'rgba(55, 128, 191, 0.7)', line: { color: 'rgba(55, 128, 191, 1.0)', width: 1 } }, opacity: 0.7 };\n", prefix, prefix, prefix);
    
    // --- Overall Fit Trace ---
    double *fit_values = malloc((hist_len + 1) * sizeof(double));
    if (!fit_values) { return; }
    calculate_theoretical_fit(fit, hist_len, fit_values);
    for (int i = 0; i < hist_len; i++) fit_values[i] *= total_counts;
    print_js_array_double(fp, "fit_y", fit_values, hist_len);
    free(fit_values);
    fprintf(fp, "var %s_fit_trace = { x: x_axis, y: fit_y, type: 'scatter', mode: 'lines', name: '%s EM Fit (%d-comp)', legendgroup: '%s', line: { color: 'red', width: 3 } };\n", prefix, view_name, fit.n_comp, prefix);

    // --- Individual Component Traces ---
    const char* component_colors[] = {"blue", "green", "orange", "purple"};
    const char* component_dash[] = {"dot", "dashdot", "longdash", "longdashdot"};
    if (fit.n_comp > 0) {
        double **comp_values = malloc(fit.n_comp * sizeof(double*));
        if (!comp_values) { return; }
        for (int i = 0; i < fit.n_comp; i++) {
            comp_values[i] = malloc((hist_len + 1) * sizeof(double));
            if (!comp_values[i]) { free(comp_values); return; }
        }
        calculate_individual_components(fit, hist_len, comp_values, total_counts);
        for (int i = 0; i < fit.n_comp; i++) {
            print_js_array_double(fp, "comp_y", comp_values[i], hist_len);
            
            const char* component_name = "";
            if (fit.p[i] == -1.0) {
                component_name = "Background (Poisson)";
            } else {
                double mean = fit.r[i] * (1.0 - fit.p[i]) / fit.p[i];
                if (i == 0 || mean < 10) component_name = "Signal (NB)";
                else component_name = "Multiplet (NB)";
            }
            
            fprintf(fp, "var %s_comp%d_trace = { x: x_axis, y: comp_y, type: 'scatter', mode: 'lines', name: 'Comp %d: %s', legendgroup: '%s', line: { color: '%s', width: 2, dash: '%s' }, opacity: 0.8 };\n", prefix, i, i+1, component_name, prefix, component_colors[i % 4], component_dash[i % 4]);
            free(comp_values[i]);
        }
        free(comp_values);
    }
    
    // --- Cutoff Traces ---
    fprintf(fp, "var %s_y_max = Math.max(...%s_hist_y) * 1.1;\n", prefix, prefix);
    fprintf(fp, "var %s_min_cutoff = { x: [%d,%d], y: [0,%s_y_max], type:'scatter', mode:'lines', name:'%s Min Signal (%d)', legendgroup: '%s', line: {color: 'darkgreen', width: 2, dash: 'dash'} };\n", prefix, fit.k_min_signal, fit.k_min_signal, prefix, view_name, fit.k_min_signal, prefix);
    fprintf(fp, "var %s_max_cutoff = { x: [%d,%d], y: [0,%s_y_max], type:'scatter', mode:'lines', name:'%s Max Signal (%d)', legendgroup: '%s', line: {color: 'darkorange', width: 2, dash: 'dash'} };\n", prefix, fit.k_max_signal, fit.k_max_signal, prefix, view_name, fit.k_max_signal, prefix);
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

static int count_features_with_reads(GArray **feature_hist, int n_features) {
    int count = 0;
    for (int i = 1; i <= n_features; i++) {
        GArray *h = feature_hist[i];
        if (h && h->len > 0) {
            long total_counts = 0;
            for (guint j = 0; j < h->len; j++) {
                total_counts += g_array_index(h, uint32_t, j);
            }
            if (total_counts > 0) {
                count++;
            }
        }
    }
    return count > 0 ? count : 1; // Avoid division by zero
}


void plot_simple_histogram(const char *directory,
                           const char *filename_leaf,
                           const char *title,
                           const char *xtitle,
                           const char *ytitle,
                           GArray *hist_data) {
    if (!hist_data || hist_data->len == 0) {
        fprintf(stderr, "Warning: Empty histogram data provided for plotting '%s'\n", title);
        return;
    }

    char filename[FILENAME_LENGTH];
    sprintf(filename, "%s/%s", directory, filename_leaf);

    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot create file %s\n", filename);
        return;
    }

    int hist_len = hist_data->len;
    uint32_t *data = (uint32_t*)hist_data->data;

    // --- Write HTML and Plotly JavaScript ---
    fprintf(fp, "<!DOCTYPE html><html><head><title>%s</title>", title);
    fprintf(fp, "<script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script></head>");
    fprintf(fp, "<body><div id=\"plot\" style=\"width:100%%;height:600px;\"></div><script>\n");

    // --- Generate Data Variables ---
    print_js_array_uint32(fp, "hist_y", data, hist_len);
    fprintf(fp, "var x_axis = [");
    for(int i=0; i<hist_len; ++i) fprintf(fp, "%d%s", i, (i==hist_len-1?"":","));
    fprintf(fp, "];\n");
    
    // --- Histogram Trace ---
    fprintf(fp, "var hist_trace = { x: x_axis, y: hist_y, type: 'bar', name: 'Frequency' };\n");

    // --- Combine Data into a single array for Plotly ---
    fprintf(fp, "var data = [hist_trace];\n");

    // --- Define Layout ---
    int last_bin_with_data = 0;
    for (int i = hist_len - 1; i >= 1; i--) {
        if (data[i] >= 1) { last_bin_with_data = i; break; }
    }
    int x_axis_max = (last_bin_with_data > 0) ? ((last_bin_with_data + 4) / 5 * 5) : 5;
     if (x_axis_max == 0) x_axis_max = 5;

    fprintf(fp, "var layout = {\n");
    fprintf(fp, "    title: { text: '%s', x: 0.5 },\n", title);
    fprintf(fp, "    xaxis: { title: '%s', type: 'linear', range: [0.5, %d] },\n", xtitle, x_axis_max);
    fprintf(fp, "    yaxis: { title: '%s', type: 'linear', autorange: true },\n", ytitle);
    fprintf(fp, "    hovermode: 'x unified',\n");
    fprintf(fp, "    updatemenus: [ \
        { buttons: [ {method: 'relayout', args: [{'yaxis.type': 'linear', 'yaxis.autorange': true}], label: 'Linear Y'}, {method: 'relayout', args: [{'yaxis.type': 'log', 'yaxis.autorange': true}], label: 'Log Y'} ], \
          direction: 'down', pad: {t: 10, r: 10}, showactive: true, type: 'dropdown', x: 0.01, xanchor: 'left', y: 1.1, yanchor: 'top' } \
    ]\n");
    fprintf(fp, "};\n");
    
    // --- Final Plotting Call ---
    fprintf(fp, "Plotly.newPlot('plot', data, layout);\n");
    fprintf(fp, "</script></body></html>\n");

    fclose(fp);
    fprintf(stderr, "Histogram plot saved to: %s\n", filename);
}

void plot_combined_histogram_with_em(const char *directory,
                                     GArray **feature_hist,
                                     NBSignalCut cumulative_em_fit,
                                     double posterior_cutoff,
                                     int n_features,
                                     double em_cumulative_limit) {

    GArray* cumulative_histogram = feature_hist[0];
    if (!cumulative_histogram || cumulative_histogram->len == 0) {
        fprintf(stderr, "Warning: Empty histogram provided for plotting\n");
        return;
    }
    if (n_features == 0) {
        fprintf(stderr, "Warning: Number of features is 0, cannot generate average histogram.\n");
        return;
    }
    
    char filename[FILENAME_LENGTH];
    sprintf(filename, "%s/umi_counts_histogram.html", directory);
    
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot create file %s\n", filename);
        return;
    }

    // --- Prepare Data for Both Plots ---
    int hist_len = cumulative_histogram->len;
    uint32_t *cumulative_hist_data = (uint32_t*)cumulative_histogram->data;

    uint32_t *avg_hist_data = malloc(hist_len * sizeof(uint32_t));
    if (!avg_hist_data) {
        fprintf(stderr, "Error: Memory allocation for average histogram failed.\n");
        fclose(fp);
        return;
    }
    long total_avg_counts = 0;
    int n_features_with_reads = count_features_with_reads(feature_hist, n_features);

    for (guint i = 0; i < hist_len; i++) {
        avg_hist_data[i] = (uint32_t)round((double)g_array_index(cumulative_histogram, uint32_t, i) / n_features_with_reads);
        total_avg_counts += avg_hist_data[i];
    }
    
    NBSignalCut avg_em_fit = cumulative_em_fit;
    avg_em_fit.total_counts_in_hist = total_avg_counts;
    determine_signal_cutoff_from_fit(&avg_em_fit, hist_len, posterior_cutoff, em_cumulative_limit);

    // --- Write HTML and Plotly JavaScript ---
    fprintf(fp, "<!DOCTYPE html><html><head><title>UMI Counts Histogram</title>");
    fprintf(fp, "<script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script></head>");
    fprintf(fp, "<body><div id=\"plot\" style=\"width:100%%;height:600px;\"></div><script>\n");

    // --- Generate Data Variables ---
    print_js_array_uint32(fp, "cumulative_hist_y", cumulative_hist_data, hist_len);
    print_js_array_uint32(fp, "avg_hist_y", avg_hist_data, hist_len);
    fprintf(fp, "var x_axis = [");
    for(int i=0; i<hist_len; ++i) fprintf(fp, "%d%s", i, (i==hist_len-1?"":","));
    fprintf(fp, "];\n");
    
    // --- Generate Trace Objects ---
    generate_traces_for_plotly(fp, "cumulative", cumulative_hist_data, hist_len, cumulative_em_fit);
    generate_traces_for_plotly(fp, "avg", avg_hist_data, hist_len, avg_em_fit);

    // --- Combine Data into a single array for Plotly ---
    fprintf(fp, "var data = [cumulative_hist_trace, cumulative_fit_trace, ");
    for(int i=0; i<cumulative_em_fit.n_comp; ++i) fprintf(fp, "cumulative_comp%d_trace, ", i);
    fprintf(fp, "cumulative_min_cutoff, cumulative_max_cutoff, ");
    
    fprintf(fp, "avg_hist_trace, avg_fit_trace, ");
    for(int i=0; i<avg_em_fit.n_comp; ++i) fprintf(fp, "avg_comp%d_trace, ", i);
    fprintf(fp, "avg_min_cutoff, avg_max_cutoff];\n");
    
    // --- Define Visibility Logic ---
    int n_cumulative_traces = 2 + cumulative_em_fit.n_comp + 2; // hist, fit, n_comp, min_cut, max_cut
    int n_avg_traces = 2 + avg_em_fit.n_comp + 2;

    // Create visibility arrays using a compatible loop
    fprintf(fp, "var cumulative_visible = [];\n");
    fprintf(fp, "for(var i=0; i<%d; ++i) cumulative_visible.push(true);\n", n_cumulative_traces);
    fprintf(fp, "for(var i=0; i<%d; ++i) cumulative_visible.push(false);\n", n_avg_traces);
    
    fprintf(fp, "var avg_visible = [];\n");
    fprintf(fp, "for(var i=0; i<%d; ++i) avg_visible.push(false);\n", n_cumulative_traces);
    fprintf(fp, "for(var i=0; i<%d; ++i) avg_visible.push(true);\n", n_avg_traces);

    // Set initial visibility to cumulative
    fprintf(fp, "for(var i=0; i<data.length; ++i) { data[i].visible = cumulative_visible[i]; }\n");

    // --- Define Layout and Controls ---
    int last_bin_with_data = 0;
    for (int i = hist_len - 1; i >= 1; i--) {
        if (cumulative_hist_data[i] >= 1) { last_bin_with_data = i; break; }
    }
    int x_axis_max = (last_bin_with_data + 4) / 5 * 5;
    if (x_axis_max == 0) x_axis_max = 5;

    fprintf(fp, "var layout = {\n");
    fprintf(fp, "    title: { text: 'Cumulative UMI Counts<br><sub>BIC: %.2f</sub>', x: 0.5 },\n", cumulative_em_fit.bic);
    fprintf(fp, "    xaxis: { title: 'UMI Counts', type: 'linear', range: [0.5, %d] },\n", x_axis_max);
    fprintf(fp, "    yaxis: { title: 'Frequency', type: 'linear', autorange: true },\n");
    fprintf(fp, "    hovermode: 'x unified',\n");
    fprintf(fp, "    legend: { traceorder: 'normal', x: 1.02, y: 1.0, xanchor: 'left', yanchor: 'top', bgcolor: 'rgba(255,255,255,0.8)', bordercolor: 'black', borderwidth: 1 },\n");
    fprintf(fp, "    updatemenus: [\n");
    fprintf(fp, "        { buttons: [ {method: 'update', args: [{'visible': cumulative_visible}, {'title.text': 'Cumulative UMI Counts<br><sub>BIC: %.2f</sub>'}], label: 'Cumulative'}, {method: 'update', args: [{'visible': avg_visible}, {'title.text': 'Average UMI Counts<br><sub>BIC: %.2f</sub>'}], label: 'Average'} ],\n", cumulative_em_fit.bic, avg_em_fit.bic);
    fprintf(fp, "          direction: 'down', pad: {t: 10, r: 10}, showactive: true, type: 'dropdown', x: 0.01, xanchor: 'left', y: 1.15, yanchor: 'top' },\n");
    fprintf(fp, "        { buttons: [ {method: 'relayout', args: [{'yaxis.type': 'linear', 'yaxis.autorange': true}], label: 'Linear Y'}, {method: 'relayout', args: [{'yaxis.type': 'log', 'yaxis.autorange': true}], label: 'Log Y'} ],\n");
    fprintf(fp, "          direction: 'down', pad: {t: 10, r: 10}, showactive: true, type: 'dropdown', x: 0.2, xanchor: 'left', y: 1.15, yanchor: 'top' }\n");
    fprintf(fp, "    ]\n};\n");
    
    // --- Final Plotting Call ---
    fprintf(fp, "Plotly.newPlot('plot', data, layout);\n");
    fprintf(fp, "</script></body></html>\n");

    fclose(fp);
    free(avg_hist_data);
    fprintf(stderr, "Combined histogram plot saved to: %s\n", filename);
}
