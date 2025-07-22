#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdint.h>


typedef enum { DIST_NB, DIST_POISSON, DIST_GAUSSIAN } DistType;

typedef double (*logf_fn)(const void *param, double x);
typedef void   (*mstep_fn)(void *param,
                           double sw, double swx, double swx2); /* ∑w, ∑wx, ∑wx² */
typedef void   (*init_fn)(void *param, double rough_mean);       /* quick seed   */
typedef double (*mean_fn)(const void *param);

typedef struct {
    size_t      param_size;    /* sizeof(param-struct)              */
    int         n_params;      /* free parameters (AIC/BIC)         */
    logf_fn     logf;          /* log-density / log-pmf             */
    mstep_fn    mstep;         /* M-step updater                    */
    init_fn     init;          /* small helper for seeding          */
    mean_fn     mean;          /* calculate mean from params        */
    const char *name;
} DistOps;
// ---------------------------------------------------------------------


// -------------------------  NB SPECIFIC  -----------------------------
typedef struct { double r, p; } NBParam;
static double nb_mean(const void *v) { const NBParam *p=v; return p->r*(1.0-p->p)/p->p; }
static double nb_logf(const void *v, double x)
{
    const NBParam *p = v;
    if (x < 0 || p->r <= 0.0 || p->p <= 0.0 || p->p >= 1.0) return -DBL_MAX;
    return lgamma(x + p->r) - lgamma(x + 1) - lgamma(p->r)
         + p->r*log(p->p) + x*log(1.0 - p->p);
}
static void nb_mstep(void *v,double sw,double swx,double swx2)
{
    NBParam *p = v;
    double mean = swx / sw;
    double var  = swx2/sw - mean*mean;
    if (var > mean) {
        p->p = mean/var;
        p->r = mean*mean/(var-mean);
        if (p->p < .01) p->p = .01;
        if (p->p > .99) p->p = .99;
        if (p->r < .1)  p->r = .1;
    }
}
static void nb_init(void *v,double m){ NBParam *p=v; p->p=.5; p->r=m; }
static const DistOps NB_OPS = { sizeof(NBParam), 2, nb_logf, nb_mstep,
                                nb_init, nb_mean, "Negative Binomial" };


// -------------------------  POISSON  ---------------------------------
typedef struct { double lambda; } PoisParam;
static double pois_mean(const void* v) { return ((PoisParam*)v)->lambda; }
static double pois_logf(const void *v,double x)
{
    const PoisParam *p = v;
    if (x < 0 || p->lambda <= 0) return -DBL_MAX;
    return x*log(p->lambda) - p->lambda - lgamma(x+1);
}
static void pois_mstep(void *v,double sw,double swx,double swx2)
{ ((PoisParam*)v)->lambda = swx / sw; }
static void pois_init(void *v,double m){ ((PoisParam*)v)->lambda = m > 0 ? m : .1; }
static const DistOps POIS_OPS = { sizeof(PoisParam), 1,
                                  pois_logf, pois_mstep, pois_init,
                                  pois_mean, "Poisson" };


// -------------------------  GAUSSIAN  --------------------------------
typedef struct { double mu, var; } GauParam;
static double gau_mean(const void* v) { return ((GauParam*)v)->mu; }
static double gau_logf(const void *v,double x)
{
    const GauParam *p=v;
    if (p->var <= 0) return -DBL_MAX;
    return -0.5*(log(2*M_PI*p->var) + (x-p->mu)*(x-p->mu)/p->var);
}
static void gau_mstep(void *v,double sw,double swx,double swx2)
{
    GauParam *p=v;
    p->mu  = swx / sw;
    p->var = swx2/sw - p->mu*p->mu;
    if (p->var < 1e-6) p->var = 1e-6;
}
static void gau_init(void *v,double m){ GauParam *p=v; p->mu=m; p->var=m ? m : 1; }
static const DistOps GAU_OPS  = { sizeof(GauParam), 2,
                                  gau_logf, gau_mstep, gau_init,
                                  gau_mean, "Gaussian" };

// ----------------------  Generic mixture  ----------------------------
typedef struct {
    /* data */
    int           n_cells, n_guides;
    int         **counts;             /* [cell][guide] */

    /* configuration */
    int           n_components;
    const DistOps **ops;              /* array of chosen distributions */

    /* parameters – flat blob ---------------------------------------- */
    void         **params_per_component; /* [comp] -> blob of size n_guides*param_size */

    double      **priors;             /* [comp][guide] */

    /* working space */
    double     ***post;               /* [cell][guide][comp] */
    double        logL;
    int            converged;
    double         bic, aic;
} MixtureEM;

/* helper to access (component,guide) parameter */
#define PARAM(em,k,g) ((void*) ((char*)(em)->params_per_component[k] + \
                       (g)*(em)->ops[k]->param_size))

/* ------------  Quick initialisation (moment-based)  --------------- */
static void seed_params(MixtureEM *em)
{
    for(int g=0; g<em->n_guides; ++g){
        /* crude guide-level mean for initialisation */
        double mean=0.0;
        int    nz  =0;
        for(int c=0;c<em->n_cells;++c){
            mean += em->counts[c][g];
            if(em->counts[c][g]) nz++;
        }
        mean = nz ? mean/nz : 0.1;

        for(int k=0;k<em->n_components;++k){
            em->ops[k]->init(PARAM(em,k,g), mean * (k+1)); /* stagger means */
            em->priors[k][g] = 1.0/em->n_components;
        }
    }
}

/* ------------  E-step (distribution-agnostic)  -------------------- */
static void e_step(MixtureEM *em)
{
    const int K = em->n_components;
    em->logL = 0.0;

    for(int i=0;i<em->n_cells;++i){
        for(int g=0;g<em->n_guides;++g){
            double max_lp = -DBL_MAX;
            double lp[K];

            for(int k=0;k<K;++k){
                lp[k] = em->ops[k]->logf(PARAM(em,k,g), em->counts[i][g])
                      + log(em->priors[k][g]);
                if(lp[k] > max_lp) max_lp = lp[k];
            }
            double sum_exp = 0.0;
            for(int k=0;k<K;++k){ sum_exp += exp(lp[k]-max_lp); }
            double log_sum = max_lp + log(sum_exp);
            em->logL += log_sum;

            for(int k=0;k<K;++k)
                em->post[i][g][k] = exp(lp[k]-log_sum);
        }
    }
}

/* ------------  M-step (delegates to plug-in)  ---------------------- */
static void m_step(MixtureEM *em)
{
    const int K = em->n_components;

    for(int g=0; g<em->n_guides; ++g){
        for(int k=0;k<K;++k){
            double sw=0, swx=0, swx2=0;
            for(int c=0;c<em->n_cells;++c){
                double w  = em->post[c][g][k];
                double x  = em->counts[c][g];
                sw   += w;
                swx  += w*x;
                swx2 += w*x*x;
            }
            em->priors[k][g] = sw / em->n_cells;
            if(sw>1e-12) em->ops[k]->mstep(PARAM(em,k,g), sw, swx, swx2);
        }
    }
}

/* ------------  Model criteria  ------------------------------------ */
static void compute_ic(MixtureEM *em)
{
    int p_per_guide = 0;
    for (int k=0; k < em->n_components; ++k)
        p_per_guide += em->ops[k]->n_params;
    p_per_guide += (em->n_components-1);          /* priors */
    int total_p     = p_per_guide * em->n_guides;
    int N           = em->n_cells * em->n_guides;

    em->aic = -2*em->logL + 2*total_p;
    em->bic = -2*em->logL + total_p*log(N);
}

/* ------------  Driver  -------------------------------------------- */
static void run_em(MixtureEM *em,int max_iter,double tol)
{
    double prev = -DBL_MAX;
    em->converged = 0;
    seed_params(em);
    for(int it=0; it<max_iter; ++it){
        e_step(em);
        if(fabs(em->logL - prev) < tol){ em->converged=1; break; }
        m_step(em);
        prev = em->logL;
    }
    compute_ic(em);
}

/* ------------  Construction / Destruction  ------------------------ */
static MixtureEM* create_em(int cells,int guides,int comp,
                            const DistOps **ops)
{
    MixtureEM *em = calloc(1,sizeof(*em));
    em->n_cells=cells; em->n_guides=guides; em->n_components=comp; em->ops=ops;

    /* parameter blob */
    em->params_per_component = malloc(comp*sizeof(void*));
    for (int k=0; k < comp; ++k) {
        em->params_per_component[k] = calloc(guides, ops[k]->param_size);
    }

    /* priors */
    em->priors = malloc(comp*sizeof(double*));
    for(int k=0;k<comp;++k){
        em->priors[k] = calloc(guides, sizeof(double));
    }

    /* posterior cube */
    em->post = malloc(cells*sizeof(double**));
    for(int i=0;i<cells;++i){
        em->post[i] = malloc(guides*sizeof(double*));
        for(int g=0;g<guides;++g)
            em->post[i][g] = calloc(comp,sizeof(double));
    }

    return em;
}

static void free_em(MixtureEM* em) {
    for (int k=0; k < em->n_components; ++k)
        free(em->params_per_component[k]);
    free(em->params_per_component);
    for (int k = 0; k < em->n_components; k++) {
        free(em->priors[k]);
    }
    free(em->priors);
    for (int i = 0; i < em->n_cells; i++) {
        for (int j = 0; j < em->n_guides; j++) {
            free(em->post[i][j]);
        }
        free(em->post[i]);
    }
    free(em->post);
    free(em);
}

/* -----------------  Public helper to fit a model  ----------------- */
MixtureEM* fit_model(int **counts,int cells,int guides,
                            int components,const DistOps **ops,
                            int max_iter,double tol)
{
    MixtureEM *em = create_em(cells,guides,components,ops);
    em->counts = counts;
    run_em(em,max_iter,tol);
    return em;
}


// -----------------------------------------------------
// PUBLIC HELPER – fit 2-vs-3 NB mixture on a *histogram*
// -----------------------------------------------------
typedef struct {
    int     n_comp;                /* chosen K (2 or 3)         */
    int     k_min_signal, k_max_signal;
    double  bic;
    /* parameters for up to 3 components                       *
     *    weight[k]  – prior π_k                               *
     *    r[k], p[k] – NB parameters                           */
    double  weight[3], r[3], p[3];
} NBSignalCut;

/* -----------------------------------------------------------------
 * Fit 1, 2, and 3-component models to a histogram of counts
 *   For 2/3 components, the first is Poisson (noise), the others are NB.
 * Decide by BIC, then calculate the contiguous range of k whose
 * posterior probability of originating from “signal” components
 * (all but the one with the smallest mean) exceeds `cut_off`.
 * ----------------------------------------------------------------*/
NBSignalCut
em_nb_signal_cut(const uint32_t *hist, int len, double cut_off,
                 int max_iter, double tol, uint16_t min_counts)
{
    /* ---------- 1. decompress hist → vector<int> counts ---------- */
    int n_cells = 0;
    for(int k=0;k<len;++k) n_cells += hist[k];
    if (n_cells == 0) {
        NBSignalCut out = {0};
        out.n_comp = 2; // Default to 2-component if no data
        out.k_min_signal = len;
        return out;
    }
    int *vec = malloc(n_cells*sizeof(int));
    int idx = 0;
    for(int k=0;k<len;++k)
        for(uint32_t f=0;f<hist[k];++f) vec[idx++] = k;

    int **cnt = malloc(n_cells*sizeof(int*));
    for(int i=0;i<n_cells;++i) cnt[i] = &vec[i];

    /* ---------- 2. fit 2 and 3-component models ----------------- */
    const DistOps* ops2[2] = { &POIS_OPS, &NB_OPS };
    const DistOps* ops3[3] = { &POIS_OPS, &NB_OPS, &NB_OPS };

    MixtureEM *em2 = fit_model(cnt,n_cells,1,2, ops2, max_iter,tol);
    MixtureEM *em3 = fit_model(cnt,n_cells,1,3, ops3, max_iter,tol);

    MixtureEM *best = (em3->bic < em2->bic) ? em3 : em2;

    NBSignalCut out = {0};
    out.bic           = best->bic;
    out.n_comp        = best->n_components;
    for(int k=0;k<out.n_comp;++k){
        out.weight[k] = best->priors[k][0];
        if (best->ops[k] == &NB_OPS) {
            NBParam *pp = (NBParam*)PARAM(best,k,0);
            out.r[k]      = pp->r;
            out.p[k]      = pp->p;
        } else {
            PoisParam *pp = (PoisParam*)PARAM(best,k,0);
            out.r[k]      = pp->lambda;
            out.p[k]      = -1.0; // Sentinel for Poisson
        }
    }

    /* ---------- 3. identify noise, signal, and multiplet components ------------- */
    int K = best->n_components;
    double means[3] = {0,0,0};
    for(int k=0;k<K;++k)
        means[k] = best->ops[k]->mean(PARAM(best,k,0));

    int bg_idx = 0;
    int multiplet_idx = -1; // Only used for K=3

    if (K == 3) {
        // In a 3-component model, identify background (lowest mean)
        // and multiplet (highest mean). The middle is signal.
        int max_idx = 0;
        for (int k = 1; k < K; k++) {
            if (means[k] < means[bg_idx]) bg_idx = k;
            if (means[k] > means[max_idx]) max_idx = k;
        }
        multiplet_idx = max_idx;
    } else { // K == 2
        // In a 2-component model, the one with the lower mean is background.
        for(int k=1;k<K;++k)
            if(means[k] < means[bg_idx]) bg_idx = k;
    }

    /* ---------- 4. posterior P(signal|k) on the grid ------------- */
    double pri[3]; for(int k=0;k<K;++k) pri[k] = best->priors[k][0];

    int first=-1,last=-1;
    for(int k=0;k<len;++k){
        double num = 0, den = 0;
        for(int j=0;j<K;++j){
            double lp = best->ops[j]->logf(PARAM(best, j, 0), k);
            double pk = exp(lp);
            den += pri[j]*pk;
            if (K == 3) {
                // Signal is the middle component (not background and not multiplet)
                if (j != bg_idx && j != multiplet_idx) num += pri[j]*pk;
            } else { // K == 2
                // Signal is the non-background component
                if (j != bg_idx) num += pri[j]*pk;
            }
        }
        double post = den? num/den : 0.0;
        if(post >= cut_off){
            if(first<0) first=k;
            last = k;
        }
    }

    out.k_min_signal  = first<0 ? len : first;
    out.k_max_signal  = last;

    /* ---------- 5. cleanup & return ------------------------------ */
    free_em(em2);
    free_em(em3);
    free(cnt);
    free(vec);
    return out;
}


/* -------------------------  Main demo  ---------------------------- */
int demo(void)
{
    const int N=400, G=10;
    /* make a tiny fake data-set */
    int **cnt = calloc(N,sizeof(int*));
    for(int i=0;i<N;++i){
        cnt[i]=calloc(G,sizeof(int));
        for(int g=0;g<G;++g){
            double r = (double)rand()/RAND_MAX;
            cnt[i][g] = r<0.8 ? rand()%4 : 10+rand()%15;
        }
    }

    /* Fit the three candidate laws */
    const DistOps* ops_nb[2] = {&NB_OPS, &NB_OPS};
    MixtureEM *em_nb  = fit_model(cnt,N,G,2,ops_nb,200,1e-6);

    const DistOps* ops_po[2] = {&POIS_OPS, &POIS_OPS};
    MixtureEM *em_po  = fit_model(cnt,N,G,2,ops_po,200,1e-6);

    const DistOps* ops_ga[2] = {&GAU_OPS, &GAU_OPS};
    MixtureEM *em_ga  = fit_model(cnt,N,G,2,ops_ga,200,1e-6);

    printf("\n%-18s  LL=%.2f  BIC=%.2f\n", NB_OPS.name , em_nb->logL,  em_nb->bic);
    printf("%-18s  LL=%.2f  BIC=%.2f\n",   POIS_OPS.name, em_po->logL, em_po->bic);
    printf("%-18s  LL=%.2f  BIC=%.2f\n",   GAU_OPS.name , em_ga->logL, em_ga->bic);

    /* … choose the best, inspect, etc … */
    free_em(em_nb);
    free_em(em_po);
    free_em(em_ga);

    for(int i=0;i<N;++i) free(cnt[i]);
    free(cnt);

    return 0;
}
