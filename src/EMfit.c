#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

// Forward declarations
typedef struct NBParams NBParams;
typedef struct MixtureEM MixtureEM;

// Negative Binomial parameters
struct NBParams {
    double r;  // dispersion
    double p;  // success probability
    double mean;  // cached mean = r(1-p)/p
};

// Complete mixture EM structure
struct MixtureEM {
    // Data
    int n_cells;
    int n_guides;
    int** counts;  // [cell][guide]

    // Model parameters
    int n_components;
    NBParams** params;      // [component][guide]
    double** priors;        // [component][guide] mixing proportions

    // Working arrays
    double*** posteriors;   // [cell][guide][component]
    double log_likelihood;

    // Model selection
    double bic;
    double aic;
    int converged;
};

// Calculate mean from NB parameters
double nb_mean(const NBParams* params) {
    return params->r * (1.0 - params->p) / params->p;
}

// Update cached mean
void nb_update_mean(NBParams* params) {
    params->mean = nb_mean(params);
}

// Robust negative binomial log PMF
double nb_logpmf(int k, const NBParams* params) {
    if (k < 0 || params->r <= 0 || params->p <= 0 || params->p >= 1) {
        return -DBL_MAX;
    }

    return lgamma(k + params->r) - lgamma(k + 1) - lgamma(params->r) +
           params->r * log(params->p) + k * log(1 - params->p);
}

// Initialize parameters using k-means style approach
void initialize_parameters(MixtureEM* em) {
    // For each guide, find min/max/median counts
    for (int g = 0; g < em->n_guides; g++) {
        double* guide_counts = (double*)malloc(em->n_cells * sizeof(double));
        int n_nonzero = 0;

        // Collect non-zero counts
        for (int c = 0; c < em->n_cells; c++) {
            if (em->counts[c][g] > 0) {
                guide_counts[n_nonzero++] = (double)em->counts[c][g];
            }
        }

        if (n_nonzero == 0) {
            // Handle guides with all zeros
            for (int k = 0; k < em->n_components; k++) {
                em->params[k][g].r = 0.1;
                em->params[k][g].p = 0.99;
                nb_update_mean(&em->params[k][g]);
                em->priors[k][g] = 1.0 / em->n_components;
            }
            free(guide_counts);
            continue;
        }

        // Sort counts
        for (int i = 0; i < n_nonzero - 1; i++) {
            for (int j = i + 1; j < n_nonzero; j++) {
                if (guide_counts[i] > guide_counts[j]) {
                    double temp = guide_counts[i];
                    guide_counts[i] = guide_counts[j];
                    guide_counts[j] = temp;
                }
            }
        }

        // Calculate quantiles for initialization
        double q25 = guide_counts[n_nonzero / 4];
        double q50 = guide_counts[n_nonzero / 2];
        double q75 = guide_counts[3 * n_nonzero / 4];

        if (em->n_components == 2) {
            // Component 0: Ambient (lower counts)
            double mean0 = q25;
            em->params[0][g].p = 0.8;
            em->params[0][g].r = mean0 * em->params[0][g].p / (1 - em->params[0][g].p);
            nb_update_mean(&em->params[0][g]);

            // Component 1: True signal (higher counts)
            double mean1 = q75;
            em->params[1][g].p = 0.3;
            em->params[1][g].r = mean1 * em->params[1][g].p / (1 - em->params[1][g].p);
            nb_update_mean(&em->params[1][g]);

            // Initial mixing proportions
            em->priors[0][g] = 0.7;
            em->priors[1][g] = 0.3;
        } else if (em->n_components == 3) {
            // Component 0: Ambient
            double mean0 = q25;
            em->params[0][g].p = 0.8;
            em->params[0][g].r = mean0 * em->params[0][g].p / (1 - em->params[0][g].p);
            nb_update_mean(&em->params[0][g]);

            // Component 1: Singlet
            double mean1 = q50;
            em->params[1][g].p = 0.4;
            em->params[1][g].r = mean1 * em->params[1][g].p / (1 - em->params[1][g].p);
            nb_update_mean(&em->params[1][g]);

            // Component 2: Doublet (approximately 2x singlet)
            double mean2 = q75;
            em->params[2][g].p = 0.3;
            em->params[2][g].r = mean2 * em->params[2][g].p / (1 - em->params[2][g].p);
            nb_update_mean(&em->params[2][g]);

            // Initial mixing proportions
            em->priors[0][g] = 0.6;
            em->priors[1][g] = 0.35;
            em->priors[2][g] = 0.05;
        }

        free(guide_counts);
    }
}

// E-step: compute posteriors
void em_e_step(MixtureEM* em) {
    em->log_likelihood = 0.0;

    for (int i = 0; i < em->n_cells; i++) {
        for (int j = 0; j < em->n_guides; j++) {
            double log_probs[em->n_components];
            double max_log_prob = -DBL_MAX;

            // Compute log probabilities for each component
            for (int k = 0; k < em->n_components; k++) {
                log_probs[k] = nb_logpmf(em->counts[i][j], &em->params[k][j]) +
                              log(em->priors[k][j]);
                if (log_probs[k] > max_log_prob) {
                    max_log_prob = log_probs[k];
                }
            }

            // Compute posteriors using log-sum-exp trick
            double sum_exp = 0.0;
            for (int k = 0; k < em->n_components; k++) {
                sum_exp += exp(log_probs[k] - max_log_prob);
            }
            double log_sum = max_log_prob + log(sum_exp);

            // Store posteriors
            for (int k = 0; k < em->n_components; k++) {
                em->posteriors[i][j][k] = exp(log_probs[k] - log_sum);
            }

            em->log_likelihood += log_sum;
        }
    }
}

// M-step: update parameters
void em_m_step(MixtureEM* em) {
    // For each guide and component
    for (int g = 0; g < em->n_guides; g++) {
        for (int k = 0; k < em->n_components; k++) {
            double sum_posterior = 0.0;
            double weighted_sum = 0.0;
            double weighted_sum_sq = 0.0;

            // Calculate sufficient statistics
            for (int c = 0; c < em->n_cells; c++) {
                double w = em->posteriors[c][g][k];
                double x = (double)em->counts[c][g];

                sum_posterior += w;
                weighted_sum += w * x;
                weighted_sum_sq += w * x * x;
            }

            // Update mixing proportion
            em->priors[k][g] = sum_posterior / em->n_cells;

            // Update NB parameters using method of moments
            if (sum_posterior > 1e-10) {
                double mean = weighted_sum / sum_posterior;
                double var = (weighted_sum_sq / sum_posterior) - mean * mean;

                // Ensure variance > mean for NB
                if (var > mean + 1e-10) {
                    em->params[k][g].p = mean / var;
                    em->params[k][g].r = mean * mean / (var - mean);

                    // Bounds checking
                    if (em->params[k][g].p <= 0) em->params[k][g].p = 0.01;
                    if (em->params[k][g].p >= 1) em->params[k][g].p = 0.99;
                    if (em->params[k][g].r <= 0) em->params[k][g].r = 0.1;
                    if (em->params[k][g].r > 1000) em->params[k][g].r = 1000;
                } else {
                    // Handle underdispersion (approximate with high r)
                    em->params[k][g].p = 0.5;
                    em->params[k][g].r = mean * 2;
                }

                nb_update_mean(&em->params[k][g]);
            }
        }
    }
}

// Calculate BIC and AIC for model selection
void calculate_model_criteria(MixtureEM* em) {
    // Number of parameters per guide: 
    // - (r, p) for each component
    // - (k-1) mixing proportions
    int params_per_guide = em->n_components * 2 + (em->n_components - 1);
    int total_params = params_per_guide * em->n_guides;
    int n_observations = em->n_cells * em->n_guides;

    em->aic = -2 * em->log_likelihood + 2 * total_params;
    em->bic = -2 * em->log_likelihood + total_params * log(n_observations);
}

// Run EM algorithm
void run_em(MixtureEM* em, int max_iter, double tol) {
    double prev_ll = -DBL_MAX;
    em->converged = 0;

    initialize_parameters(em);

    for (int iter = 0; iter < max_iter; iter++) {
        // E-step
        em_e_step(em);

        // Check convergence
        double ll_change = fabs(em->log_likelihood - prev_ll);
        if (ll_change < tol) {
            em->converged = 1;
            printf("Converged after %d iterations (LL change: %.6f)\n", 
                   iter + 1, ll_change);
            break;
        }

        // M-step
        em_m_step(em);

        prev_ll = em->log_likelihood;

        if ((iter + 1) % 10 == 0) {
            printf("Iteration %d: LL = %.2f\n", iter + 1, em->log_likelihood);
        }
    }

    calculate_model_criteria(em);
}

// Create EM structure
MixtureEM* create_mixture_em(int n_cells, int n_guides, int n_components) {
    MixtureEM* em = (MixtureEM*)calloc(1, sizeof(MixtureEM));

    em->n_cells = n_cells;
    em->n_guides = n_guides;
    em->n_components = n_components;

    // Allocate parameters
    em->params = (NBParams**)calloc(n_components, sizeof(NBParams*));
    em->priors = (double**)calloc(n_components, sizeof(double*));

    for (int k = 0; k < n_components; k++) {
        em->params[k] = (NBParams*)calloc(n_guides, sizeof(NBParams));
        em->priors[k] = (double*)calloc(n_guides, sizeof(double));
    }

    // Allocate posteriors
    em->posteriors = (double***)calloc(n_cells, sizeof(double**));
    for (int i = 0; i < n_cells; i++) {
        em->posteriors[i] = (double**)calloc(n_guides, sizeof(double*));
        for (int j = 0; j < n_guides; j++) {
            em->posteriors[i][j] = (double*)calloc(n_components, sizeof(double));
        }
    }

    return em;
}

// Model selection: fit both 2 and 3 component models and choose best
MixtureEM* fit_best_model(int** counts, int n_cells, int n_guides,
                         int max_iter, double tol) {
    printf("Fitting 2-component model...\n");
    MixtureEM* em2 = create_mixture_em(n_cells, n_guides, 2);
    em2->counts = counts;
    run_em(em2, max_iter, tol);

    printf("\nFitting 3-component model...\n");
    MixtureEM* em3 = create_mixture_em(n_cells, n_guides, 3);
    em3->counts = counts;
    run_em(em3, max_iter, tol);

    printf("\nModel Selection:\n");
    printf("2-component: LL=%.2f, AIC=%.2f, BIC=%.2f\n", 
           em2->log_likelihood, em2->aic, em2->bic);
    printf("3-component: LL=%.2f, AIC=%.2f, BIC=%.2f\n", 
           em3->log_likelihood, em3->aic, em3->bic);

    // Choose model with lower BIC
    if (em2->bic < em3->bic) {
        printf("Selected: 2-component model (lower BIC)\n");
        // Free em3
        free_mixture_em(em3);
        return em2;
    } else {
        printf("Selected: 3-component model (lower BIC)\n");
        // Free em2
        free_mixture_em(em2);
        return em3;
    }
}

// Get guide assignments for a cell
void classify_cell_guides(MixtureEM* em, int cell_idx, double threshold) {
    printf("\nCell %d guide assignments:\n", cell_idx);

    for (int g = 0; g < em->n_guides; g++) {
        // Find component with highest posterior
        int best_comp = 0;
        double max_post = 0.0;

        for (int k = 0; k < em->n_components; k++) {
            if (em->posteriors[cell_idx][g][k] > max_post) {
                max_post = em->posteriors[cell_idx][g][k];
                best_comp = k;
            }
        }

        // Only report if count > 0 and posterior > threshold
        if (em->counts[cell_idx][g] > 0 && max_post > threshold) {
            const char* comp_names[] = {"Ambient", "Signal/Singlet", "Doublet"};
            printf("  Guide %d: count=%d, assigned to %s (posterior=%.3f)\n",
                   g, em->counts[cell_idx][g], 
                   comp_names[best_comp], max_post);
        }
    }
}

// Summary statistics
void print_em_summary(MixtureEM* em) {
    printf("\nFitted Model Summary:\n");
    printf("Components: %d\n", em->n_components);
    printf("Log-likelihood: %.2f\n", em->log_likelihood);
    printf("BIC: %.2f\n", em->bic);
    printf("AIC: %.2f\n\n", em->aic);

    // Average parameters across guides
    for (int k = 0; k < em->n_components; k++) {
        double avg_mean = 0.0;
        double avg_disp = 0.0;
        double avg_prior = 0.0;

        for (int g = 0; g < em->n_guides; g++) {
            avg_mean += em->params[k][g].mean;
            avg_disp += em->params[k][g].r;
            avg_prior += em->priors[k][g];
        }

        avg_mean /= em->n_guides;
        avg_disp /= em->n_guides;
        avg_prior /= em->n_guides;

        printf("Component %d:\n", k);
        printf("  Average mean: %.2f\n", avg_mean);
        printf("  Average dispersion: %.2f\n", avg_disp);
        printf("  Average prior: %.3f\n", avg_prior);
    }
}

// Free memory
void free_mixture_em(MixtureEM* em) {
    for (int k = 0; k < em->n_components; k++) {
        free(em->params[k]);
        free(em->priors[k]);
    }
    free(em->params);
    free(em->priors);

    for (int i = 0; i < em->n_cells; i++) {
        for (int j = 0; j < em->n_guides; j++) {
            free(em->posteriors[i][j]);
        }
        free(em->posteriors[i]);
    }
    free(em->posteriors);

    free(em);
}

// Example usage
int main() {
    int n_cells = 1000;
    int n_guides = 20;

    // Allocate and populate count matrix (example with simulated data)
    int** counts = (int**)calloc(n_cells, sizeof(int*));
    for (int i = 0; i < n_cells; i++) {
        counts[i] = (int*)calloc(n_guides, sizeof(int));

        // Simulate some data
        for (int j = 0; j < n_guides; j++) {
            double r = (double)rand() / RAND_MAX;
            if (r < 0.7) {
                counts[i][j] = rand() % 3;  // Ambient
            } else if (r < 0.95) {
                counts[i][j] = 5 + rand() % 20;  // Signal
            } else {
                counts[i][j] = 20 + rand() % 40;  // Doublet
            }
        }
    }

    // Fit best model
    MixtureEM* best_model = fit_best_model(counts, n_cells, n_guides, 100, 1e-6);

    // Print summary
    print_em_summary(best_model);

    // Classify some cells
    for (int i = 0; i < 5; i++) {
        classify_cell_guides(best_model, i, 0.5);
    }

    // Cleanup
    free_mixture_em(best_model);
    for (int i = 0; i < n_cells; i++) {
        free(counts[i]);
    }
    free(counts);

    return 0;
}
