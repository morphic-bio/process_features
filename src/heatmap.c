#include <cairo.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <hdf5.h>
#include "plasma_colormap_16.h"
#include "plasma_colormap_64.h"
#include "plasma_colormap_256.h"
#include "plasma_colormap_1024.h"
#define CELL_SIZE 10  // Size of each cell in the heatmap
#define BAR_WIDTH 20  // Width of the color bar
#define BASE_PADDING 10  // Base padding for additional adjustments
#define BAR_GRAPH_HEIGHT 100   // Height of the bar graph area
typedef struct {
    int n_obs;               // Number of rows
    int n_vars;              // Number of columns
    double **matrix;         // 2D array to hold the matrix values
    char **obs_names;        // Array to hold row (obs) names
    char **var_names;        // Array to hold column (var) names
} H5ADData;
void check_hdf5_error(hid_t status, const char *message) {
    if (status < 0) {
        fprintf(stderr, "HDF5 Error: %s\n", message);
        exit(EXIT_FAILURE);
    }
}

char **read_hdf5_string_array(hid_t file_id, const char *dataset_name, int *num_elements) {
    hid_t dataset_id = H5Dopen(file_id, dataset_name, H5P_DEFAULT);
    check_hdf5_error(dataset_id, "Failed to open dataset");

    hid_t dataspace_id = H5Dget_space(dataset_id);
    hsize_t dims[1];
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
    *num_elements = (int)dims[0];

    // Read the string data
    char **string_array = (char **)malloc(*num_elements * sizeof(char *));
    for (int i = 0; i < *num_elements; i++) {
        string_array[i] = (char *)malloc(1024 * sizeof(char)); // Allocate memory for each string
    }

    hid_t memtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(memtype, H5T_VARIABLE);
    H5Dread(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, string_array);

    // Close HDF5 resources
    H5Tclose(memtype);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);

    return string_array;
}

double **read_hdf5_matrix(hid_t file_id, const char *dataset_name, int *n_obs, int *n_vars) {
    hid_t dataset_id = H5Dopen(file_id, dataset_name, H5P_DEFAULT);
    check_hdf5_error(dataset_id, "Failed to open dataset");

    // Get the dimensions of the matrix
    hid_t dataspace_id = H5Dget_space(dataset_id);
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
    *n_obs = (int)dims[0];
    *n_vars = (int)dims[1];

    // Allocate memory for the matrix
    double **matrix = (double **)malloc(*n_obs * sizeof(double *));
    double *data = (double *)malloc((*n_obs) * (*n_vars) * sizeof(double));

    for (int i = 0; i < *n_obs; i++) {
        matrix[i] = data + i * (*n_vars);
    }

    // Read the matrix data
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    // Close HDF5 resources
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);

    return matrix;
}

H5ADData *read_h5ad(const char *filename) {
    H5ADData *h5ad_data = (H5ADData *)malloc(sizeof(H5ADData));

    // Open the file
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    check_hdf5_error(file_id, "Failed to open file");

    // Read row (obs) names
    h5ad_data->obs_names = read_hdf5_string_array(file_id, "obs/_index", &(h5ad_data->n_obs));

    // Read column (var) names
    h5ad_data->var_names = read_hdf5_string_array(file_id, "var/_index", &(h5ad_data->n_vars));

    // Read the matrix
    h5ad_data->matrix = read_hdf5_matrix(file_id, "X/data", &(h5ad_data->n_obs), &(h5ad_data->n_vars));

    // Close the HDF5 file
    H5Fclose(file_id);

    return h5ad_data;
}

void free_h5ad_data(H5ADData *h5ad_data) {
    for (int i = 0; i < h5ad_data->n_obs; i++) {
        free(h5ad_data->obs_names[i]);
    }
    for (int i = 0; i < h5ad_data->n_vars; i++) {
        free(h5ad_data->var_names[i]);
    }
    free(h5ad_data->obs_names);
    free(h5ad_data->var_names);
    free(h5ad_data->matrix[0]); // Free the flattened matrix
    free(h5ad_data->matrix);
    free(h5ad_data);
}

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
    // Clamp the normalized value to [0, 1]
    if (normalized_value < 0.0) normalized_value = 0.0;
    if (normalized_value > 1.0) normalized_value = 1.0;

    // Scale the normalized value to the colormap size
    int index = (int)(normalized_value * (colormap_size - 1)); // Map to 0–(colormap_size-1)

    // Map to RGB values from the provided colormap
    *r = colormap[index][0];
    *g = colormap[index][1];
    *b = colormap[index][2];
}


// Function to normalize values for coloring
double normalize(int value, int max_value) {
    return (double)value / max_value;
}

/* // Function to generate heatmap with dynamic padding
void generate_heatmap(const char *directory, feature_arrays *features, int **coexpression_histograms) {
    char output_file[1024];
    if (directory[strlen(directory) - 1] == '/') {
        snprintf(output_file, sizeof(output_file), "%sheatmap.png", directory);
    } else {
        snprintf(output_file, sizeof(output_file), "%s/heatmap.png", directory);
    }
    //find max_name_length
    int max_name_length=0;
    const int num_rows = features->number_of_features;
    int column_sums[num_rows];
    memset(column_sums,0,num_rows*sizeof(int));    
    for (int i=0; i<num_rows; i++){
        if (strlen(features->feature_names[i])>max_name_length){
            max_name_length=strlen(features->feature_names[i]);
        }
    }

    // Calculate dynamic padding
    int left_padding = BASE_PADDING + max_name_length * 7; // 7 pixels per character as a rough estimate
    int right_padding = BASE_PADDING + 50; // Allow space for color bar labels
    int filtered_rows = 0;
    int filter_mask[num_rows];
    int num_cols = 0;
    
    //create a filter mask - necessary later for color map and calculate max column
    //remember that coexpression_histograms is 1 indexed with the 0 column saved for totals
    for (int i = 0; i < num_rows; i++) {
       filter_mask[i]=coexpression_histograms[i+1][0]>0;
         if (filter_mask[i]){
              filtered_rows++;
         }
       if (coexpression_histograms[i+1][0]> num_cols){
           num_cols=coexpression_histograms[i+1][0];
       }
    }
    // calculate column sums and find the maximum column sum
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            column_sums[j] += coexpression_histograms[i+1][j+1];
        }
    }
    //correct for multiple counting of column sums
    for (int j = 0; j < num_cols; j++) {
        column_sums[j] /= j+1;
    }
    int max_column_sum = 0;
    for (int j = 0; j < num_cols; j++) {
        if (column_sums[j] > max_column_sum) {
            max_column_sum = column_sums[j];
        }
    }
    // Adjust canvas dimensions to include the bar graph
    int width = left_padding + num_cols * CELL_SIZE + BAR_WIDTH + right_padding;
    int height = BASE_PADDING + BAR_GRAPH_HEIGHT + BASE_PADDING + 20 + filtered_rows * CELL_SIZE + BASE_PADDING; // +20 for x-axis labels

    cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
    cairo_t *cr = cairo_create(surface);

    // Draw bar graph
    for (int j = 0; j < num_cols; j++) {
        double bar_height = (double)column_sums[j] / max_column_sum * BAR_GRAPH_HEIGHT;
        cairo_set_source_rgb(cr, 0.2, 0.4, 0.8); // Bar color (blue)
        cairo_rectangle(cr, left_padding + j * CELL_SIZE, BASE_PADDING + BAR_GRAPH_HEIGHT - bar_height, CELL_SIZE, bar_height);
        cairo_fill(cr);
    }
    // **Draw side labels for maximum and midpoint**
    cairo_set_font_size(cr, 10);
    cairo_set_source_rgb(cr, 0, 0, 0); // Black color for labels

    // Maximum value label
    char max_label[20];
    snprintf(max_label, sizeof(max_label), "%-d", max_column_sum);
    cairo_move_to(cr, left_padding - 30 - 10, BASE_PADDING + 5); // Position label at top-left of bar graph
    cairo_show_text(cr, max_label);

    // Midpoint value label
    char midpoint_label[20];
    int midpoint_value = max_column_sum / 2;
    snprintf(midpoint_label, sizeof(midpoint_label), "%d", midpoint_value);
    cairo_move_to(cr, left_padding - 30 - 10, BASE_PADDING + BAR_GRAPH_HEIGHT / 2 + 5); // Position at the midpoint of the bar graph
    cairo_show_text(cr, midpoint_label);
    
    // Draw x-axis labels under the bar graph
    cairo_set_font_size(cr, 10);
    cairo_set_source_rgb(cr, 0, 0, 0); // Black color for labels
    for (int j = 0; j < num_cols; j++) {
        if ((j+1) % 5 == 0) { // Label every 5th column
            char label[10];
            snprintf(label, sizeof(label), "%d", j+1);
            cairo_move_to(cr, left_padding + j * CELL_SIZE + CELL_SIZE / 4, BASE_PADDING + BAR_GRAPH_HEIGHT + 15); // Adjust placement
            cairo_show_text(cr, label);
        }
    }

    // Draw heatmap
    //find max value
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
            //fprintf(stderr, "Intensity %f\n", intensity);
            double r, g, b;
            value_to_color(intensity, &r, &g, &b, plasma_colormap, colormap_size);
            cairo_set_source_rgb(cr, r, g, b);
            cairo_rectangle(cr, left_padding + j * CELL_SIZE, BASE_PADDING + BAR_GRAPH_HEIGHT + 20 + draw_row * CELL_SIZE, CELL_SIZE, CELL_SIZE);
            cairo_fill(cr);
        }
        // Draw row labels
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 10);
        cairo_move_to(cr, BASE_PADDING, BASE_PADDING + BAR_GRAPH_HEIGHT + 20 + draw_row * CELL_SIZE + CELL_SIZE / 2);
        cairo_show_text(cr, features->feature_names[i]);
        draw_row++;
    }
    // Draw color bar
    for (int i = 0; i < height - BAR_GRAPH_HEIGHT - 20 - BASE_PADDING*3 ; i++) {
        double normalized_value = 1.0 - (double)i / (height - BAR_GRAPH_HEIGHT - BASE_PADDING * 3 - 20);
        double r, g, b;
        value_to_color(normalized_value, &r, &g, &b, plasma_colormap, colormap_size);
        cairo_set_source_rgb(cr, r, g, b);
        cairo_rectangle(cr, width - right_padding - BAR_WIDTH, BASE_PADDING + BAR_GRAPH_HEIGHT + 20 + i, BAR_WIDTH, 1);
        cairo_fill(cr);
    }

    // Add labels to the color bar
    cairo_set_font_size(cr, 10);
    for (int i = 0; i <= height - BAR_GRAPH_HEIGHT - BASE_PADDING * 3 - 20; i += 10 * CELL_SIZE) {
        double normalized_value = 1.0 - (double)i / (height - BAR_GRAPH_HEIGHT - BASE_PADDING * 3 - 20 );
        int value = (int)(normalized_value * max_value);
        char label[20];
        snprintf(label, sizeof(label), "%d", value);
        cairo_move_to(cr, width - right_padding + 5, BASE_PADDING + BAR_GRAPH_HEIGHT + 20 + i + 4); // Adjust for padding
        cairo_show_text(cr, label);
    }

    // Save to PNG
    cairo_surface_write_to_png(surface, output_file);

    // Clean up
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
} */


