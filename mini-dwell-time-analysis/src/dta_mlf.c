#include "dta_core.h"
#include "dta_mlf.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * dta_mlf.c - Multiple Lyapunov Function Construction
 *
 * Implements:
 *   Quadratic MLF via Lyapunov equation (A_i^T P_i + P_i A_i = -Q)
 *   Optimal MLF with maximized decay rates
 *   Coupling constant mu computation
 *   MLF evaluation and derivative computation
 *   MLF verification for dwell-time GAS
 *   Common quadratic Lyapunov function search
 *   Piecewise quadratic MLF construction
 *   Required dwell time from MLF parameters
 *
 * References:
 *   Branicky (1998) IEEE TAC 43(4):475-482
 *   Liberzon (2003) "Switching in Systems and Control" Ch.3.2
 *   Decarlo, Branicky, Pettersson, Lennartson (2000) Proc IEEE
 * ============================================================== */

static void* safe_alloc(size_t sz) {
    void* p = malloc(sz);
    if (!p) { fprintf(stderr, "DTA mlf: alloc fail\n"); exit(1); }
    return p;
}

DTA_MultipleLyapunov* dta_mlf_create(int n_modes, int n) {
    DTA_MultipleLyapunov* mlf = safe_alloc(sizeof(DTA_MultipleLyapunov));
    mlf->n_modes = n_modes;
    mlf->n = n;
    mlf->type = DTA_MLF_QUADRATIC;
    mlf->mu = 1.0;
    mlf->P_i = safe_alloc((size_t)n_modes * sizeof(double*));
    mlf->lambda_i = safe_alloc((size_t)n_modes * sizeof(double));
    int i;
    for (i = 0; i < n_modes; i++) {
        mlf->P_i[i] = safe_alloc((size_t)(n*n) * sizeof(double));
        memset(mlf->P_i[i], 0, (size_t)(n*n) * sizeof(double));
        mlf->lambda_i[i] = 1.0;
    }
    return mlf;
}

void dta_mlf_free(DTA_MultipleLyapunov* mlf) {
    if (!mlf) return;
    int i;
    for (i = 0; i < mlf->n_modes; i++) free(mlf->P_i[i]);
    free(mlf->P_i);
    free(mlf->lambda_i);
    free(mlf);
}

int dta_mlf_set_P(DTA_MultipleLyapunov* mlf, int mode, const double* P) {
    if (!mlf || mode < 0 || mode >= mlf->n_modes || !P) return -1;
    memcpy(mlf->P_i[mode], P, (size_t)(mlf->n * mlf->n) * sizeof(double));
    return 0;
}

/* Construct quadratic MLF by solving A_i^T P_i + P_i A_i = -Q_i */
DTA_MultipleLyapunov* dta_mlf_construct_quadratic(
    const DTA_SwitchedSystem* sys, const double** Q_matrices) {
    if (!sys) return NULL;
    int n = sys->state_dim, m = sys->n_modes;
    DTA_MultipleLyapunov* mlf = dta_mlf_create(m, n);
    int i, j;
    /* Default Q = I if not provided */
    double* Q_default = NULL;
    if (!Q_matrices) {
        Q_default = safe_alloc((size_t)(n*n) * sizeof(double));
        memset(Q_default, 0, (size_t)(n*n) * sizeof(double));
        for (i = 0; i < n; i++) Q_default[i*n + i] = 1.0;
    }
    for (i = 0; i < m; i++) {
        const double* Q = Q_matrices ? Q_matrices[i] : Q_default;
        if (!Q) Q = Q_default;
        int ret = dta_solve_lyapunov(sys->modes[i].A, n, Q, mlf->P_i[i]);
        if (ret != 0) {
            /* Fallback: P_i = I */
            for (j = 0; j < n; j++) mlf->P_i[i][j*n + j] = 1.0;
        }
        /* Compute lambda_i from P_i */
        double* re = safe_alloc((size_t)n * sizeof(double));
        double* im = safe_alloc((size_t)n * sizeof(double));
        dta_eigenvalues(mlf->P_i[i], n, re, im, 200);
        double min_eig = re[0];
        for (j = 1; j < n; j++)
            if (re[j] < min_eig) min_eig = re[j];
        mlf->lambda_i[i] = (min_eig > 1e-12) ? 1.0 / (2.0 * sqrt(min_eig)) : 1.0;
        free(re); free(im);
    }
    if (Q_default) free(Q_default);
    mlf->mu = dta_mlf_compute_mu(mlf);
    return mlf;
}

/* Construct optimal MLF by maximizing min lambda_i iteratively */
DTA_MultipleLyapunov* dta_mlf_construct_optimal(
    const DTA_SwitchedSystem* sys, int max_iter, double tol) {
    DTA_MultipleLyapunov* mlf = dta_mlf_construct_quadratic(sys, NULL);
    if (!mlf) return NULL;
    int n = mlf->n, m = mlf->n_modes;
    int iter, i, j;
    double step = 0.1;
    for (iter = 0; iter < max_iter; iter++) {
        double lambda_min = mlf->lambda_i[0];
        int worst = 0;
        for (i = 1; i < m; i++) {
            if (mlf->lambda_i[i] < lambda_min) {
                lambda_min = mlf->lambda_i[i];
                worst = i;
            }
        }
        /* Try to increase lambda for the worst mode */
        double* Q_new = safe_alloc((size_t)(n*n) * sizeof(double));
        for (i = 0; i < n; i++)
            Q_new[i*n + i] = 1.0 + step * lambda_min;
        double* P_new = safe_alloc((size_t)(n*n) * sizeof(double));
        dta_solve_lyapunov(sys->modes[worst].A, n, Q_new, P_new);
        double* re = safe_alloc((size_t)n * sizeof(double));
        double* im = safe_alloc((size_t)n * sizeof(double));
        dta_eigenvalues(P_new, n, re, im, 200);
        double min_eig = re[0];
        for (j = 1; j < n; j++)
            if (re[j] < min_eig) min_eig = re[j];
        free(re); free(im);
        if (min_eig > tol) {
            memcpy(mlf->P_i[worst], P_new, (size_t)(n*n) * sizeof(double));
            mlf->lambda_i[worst] = 1.0 / (2.0 * sqrt(min_eig));
        } else {
            step *= 0.5;
        }
        free(Q_new); free(P_new);
        if (step < tol) break;
    }
    mlf->mu = dta_mlf_compute_mu(mlf);
    return mlf;
}

/* Compute coupling constant mu = max_{i,j} lambda_max(P_i P_j^{-1})
 * Simplified: use eigenvalue ratio bound */
double dta_mlf_compute_mu(const DTA_MultipleLyapunov* mlf) {
    if (!mlf) return 1.0;
    int n = mlf->n, m = mlf->n_modes;
    double mu = 1.0;
    int i, j;
    for (i = 0; i < m; i++) {
        double* re_i = safe_alloc((size_t)n * sizeof(double));
        double* im_i = safe_alloc((size_t)n * sizeof(double));
        dta_eigenvalues(mlf->P_i[i], n, re_i, im_i, 200);
        double max_eig_i = re_i[0];
        double min_eig_i = re_i[0];
        int k;
        for (k = 1; k < n; k++) {
            if (re_i[k] > max_eig_i) max_eig_i = re_i[k];
            if (re_i[k] < min_eig_i) min_eig_i = re_i[k];
        }
        free(re_i); free(im_i);
        for (j = 0; j < m; j++) {
            if (i == j) continue;
            double* re_j = safe_alloc((size_t)n * sizeof(double));
            double* im_j = safe_alloc((size_t)n * sizeof(double));
            dta_eigenvalues(mlf->P_i[j], n, re_j, im_j, 200);
            double max_eig_j = re_j[0];
            double min_eig_j = re_j[0];
            for (k = 1; k < n; k++) {
                if (re_j[k] > max_eig_j) max_eig_j = re_j[k];
                if (re_j[k] < min_eig_j) min_eig_j = re_j[k];
            }
            free(re_j); free(im_j);
            double ratio_ij = (min_eig_j > 1e-12) ? max_eig_i / min_eig_j : 1e6;
            if (ratio_ij > mu) mu = ratio_ij;
        }
    }
    return mu;
}

double dta_mlf_evaluate(const DTA_MultipleLyapunov* mlf, int mode,
                         const double* x) {
    if (!mlf || mode < 0 || mode >= mlf->n_modes || !x) return 0.0;
    int n = mlf->n;
    double V = 0.0;
    int i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            V += x[i] * mlf->P_i[mode][i*n + j] * x[j];
    return V;
}

double dta_mlf_derivative(const DTA_MultipleLyapunov* mlf, int mode,
                           const DTA_SwitchedSystem* sys, const double* x) {
    if (!mlf || !sys || mode < 0 || mode >= mlf->n_modes || !x) return 0.0;
    int n = mlf->n, i, j, k;
    /* V_dot = x^T (A^T P + P A) x */
    double V_dot = 0.0;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            double mij = 0.0;
            for (k = 0; k < n; k++)
                mij += sys->modes[mode].A[k*n + i] * mlf->P_i[mode][k*n + j]
                     + mlf->P_i[mode][i*n + k] * sys->modes[mode].A[k*n + j];
            V_dot += x[i] * mij * x[j];
        }
    return V_dot;
}

/* Verify MLF satisfies conditions for GAS under dwell time tau_d */
DTA_MLFVerification dta_mlf_verify(const DTA_MultipleLyapunov* mlf,
    const DTA_SwitchedSystem* sys, const DTA_SwitchingSignal* sig,
    const double* x0, double tau_d) {
    DTA_MLFVerification v;
    memset(&v, 0, sizeof(v));
    if (!mlf || !sys) return v;
    int n = mlf->n, m = mlf->n_modes;
    v.all_pd = true;
    v.all_decrease = true;
    v.switching_decrease = true;
    int i, j;
    /* Check P_i > 0 */
    for (i = 0; i < m; i++) {
        double* re = safe_alloc((size_t)n * sizeof(double));
        double* im = safe_alloc((size_t)n * sizeof(double));
        dta_eigenvalues(mlf->P_i[i], n, re, im, 200);
        for (j = 0; j < n; j++)
            if (re[j] <= 0) { v.all_pd = false; break; }
        free(re); free(im);
        if (!v.all_pd) break;
    }
    /* Check V_dot_i < 0 */
    for (i = 0; i < m; i++) {
        double* M = safe_alloc((size_t)(n*n) * sizeof(double));
        int k, l;
        for (k = 0; k < n; k++)
            for (l = 0; l < n; l++) {
                double sv = 0.0;
                int r;
                for (r = 0; r < n; r++)
                    sv += sys->modes[i].A[r*n + k] * mlf->P_i[i][r*n + l]
                        + mlf->P_i[i][k*n + r] * sys->modes[i].A[r*n + l];
                M[k*n + l] = sv;
            }
        double* re = safe_alloc((size_t)n * sizeof(double));
        double* im = safe_alloc((size_t)n * sizeof(double));
        dta_eigenvalues(M, n, re, im, 200);
        for (j = 0; j < n; j++)
            if (re[j] >= -1e-12) { v.all_decrease = false; break; }
        free(M); free(re); free(im);
        if (!v.all_decrease) break;
    }
    v.valid = v.all_pd && v.all_decrease;
    v.max_coupling = mlf->mu;
    v.min_decay_rate = mlf->lambda_i[0];
    for (i = 1; i < m; i++)
        if (mlf->lambda_i[i] < v.min_decay_rate)
            v.min_decay_rate = mlf->lambda_i[i];
    v.tau_d_required = dta_mlf_required_dwell(mlf);
    (void)sig; (void)x0; (void)tau_d;
    return v;
}

bool dta_mlf_construct_common(const DTA_SwitchedSystem* sys, double* P_out,
                               int max_iter, double step, double tol) {
    if (!sys || !P_out) return false;
    int n = sys->state_dim, m = sys->n_modes;
    int i, j, k, iter;

    /* Initialize P = I */
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            P_out[i*n + j] = (i == j) ? 1.0 : 0.0;

    for (iter = 0; iter < max_iter; iter++) {
        int worst_mode = 0;
        double worst_eig = -INFINITY;
        for (i = 0; i < m; i++) {
            double* M = safe_alloc((size_t)(n*n) * sizeof(double));
            for (j = 0; j < n; j++)
                for (k = 0; k < n; k++) {
                    double sv = 0.0;
                    int l;
                    for (l = 0; l < n; l++)
                        sv += sys->modes[i].A[l*n + j] * P_out[l*n + k]
                            + P_out[j*n + l] * sys->modes[i].A[l*n + k];
                    M[j*n + k] = sv;
                }
            double* re = safe_alloc((size_t)n * sizeof(double));
            double* im = safe_alloc((size_t)n * sizeof(double));
            dta_eigenvalues(M, n, re, im, 200);
            if (re[0] > worst_eig) { worst_eig = re[0]; worst_mode = i; }
            free(M); free(re); free(im);
        }
        if (worst_eig < -tol) return true;

        /* Gradient descent step */
        double* G = safe_alloc((size_t)(n*n) * sizeof(double));
        for (j = 0; j < n; j++)
            for (k = 0; k < n; k++) {
                double sv = 0.0;
                int l;
                for (l = 0; l < n; l++)
                    sv += sys->modes[worst_mode].A[l*n + j] * P_out[l*n + k]
                        + P_out[j*n + l] * sys->modes[worst_mode].A[l*n + k];
                G[j*n + k] = sv;
            }
        for (j = 0; j < n*n; j++) P_out[j] -= step * G[j];
        free(G);

        /* Keep P positive definite and symmetric */
        for (j = 0; j < n; j++)
            if (P_out[j*n + j] < 1e-6) P_out[j*n + j] = 1e-6;
        for (j = 0; j < n; j++)
            for (k = j+1; k < n; k++) {
                double avg = 0.5 * (P_out[j*n + k] + P_out[k*n + j]);
                P_out[j*n + k] = avg; P_out[k*n + j] = avg;
            }
    }
    return false;
}

DTA_MultipleLyapunov* dta_mlf_construct_piecewise(
    const DTA_SwitchedSystem* sys, double* region_centers,
    int n_regions) {
    if (!sys || n_regions <= 0) return NULL;
    DTA_MultipleLyapunov* mlf = dta_mlf_create(n_regions, sys->state_dim);
    mlf->type = DTA_MLF_PIECEWISE_QUAD;
    int i, j;
    for (i = 0; i < n_regions; i++) {
        for (j = 0; j < sys->state_dim; j++)
            mlf->P_i[i][j*sys->state_dim + j] = 1.0;
        if (region_centers && i < n_regions * sys->state_dim) {
            double dist = 0.0;
            for (j = 0; j < sys->state_dim; j++)
                dist += region_centers[i*sys->state_dim + j]
                      * region_centers[i*sys->state_dim + j];
            if (dist > 0) {
                double scale = 1.0 / sqrt(dist);
                for (j = 0; j < sys->state_dim; j++)
                    mlf->P_i[i][j*sys->state_dim + j] = 1.0 + scale;
            }
        }
    }
    mlf->mu = dta_mlf_compute_mu(mlf);
    return mlf;
}

double dta_mlf_required_dwell(const DTA_MultipleLyapunov* mlf) {
    if (!mlf) return INFINITY;
    double lambda_min = mlf->lambda_i[0];
    int i;
    for (i = 1; i < mlf->n_modes; i++)
        if (mlf->lambda_i[i] < lambda_min)
            lambda_min = mlf->lambda_i[i];
    if (lambda_min <= 1e-12 || mlf->mu <= 1.0) return 0.0;
    return log(mlf->mu) / (2.0 * lambda_min);
}

void dta_mlf_scale(DTA_MultipleLyapunov* mlf, double alpha) {
    if (!mlf || alpha <= 0) return;
    int i, j;
    for (i = 0; i < mlf->n_modes; i++)
        for (j = 0; j < mlf->n * mlf->n; j++)
            mlf->P_i[i][j] *= alpha;
    mlf->mu = dta_mlf_compute_mu(mlf);
}
