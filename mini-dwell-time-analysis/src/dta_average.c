#include "dta_core.h"
#include "dta_switch_signal.h"
#include "dta_average.h"
#include "dta_stability.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * dta_average.c - Average Dwell-Time Theory
 *
 * Implements:
 *   Minimum ADT computation from Lyapunov data
 *   Mode-dependent ADT computation
 *   ADT signal validation (N_sigma bound check)
 *   ADT stability analysis
 *   Chatter bound N0 computation from signal
 *   ADT-constrained signal generation
 *   Convergence rate bounds under ADT
 *   Constant vs average dwell time comparison
 *
 * Theorem (Hespanha & Morse 1999):
 *   If tau_a > tau_a* = ln(c)/lambda, then GUES.
 *   If mode-dependent: tau_{a,i} > ln(mu_i)/lambda_i.
 *
 * References:
 *   Hespanha & Morse (1999) IEEE CDC
 *   Hespanha (2004) IEEE TAC 49(4):470-482
 *   Zhao & Hill (2008) Automatica 44(7):1809-1816
 * ============================================================== */

static void* safe_alloc(size_t sz) {
    void* p = malloc(sz);
    if (!p) { fprintf(stderr, "DTA adt: alloc fail\n"); exit(1); }
    return p;
}

/* Compute minimum ADT from system data.
 * tau_a* = max_i { ln(cond(P_i)^{1/2}) / lambda_i }
 * where cond(P_i) = lambda_max(P_i) / lambda_min(P_i) */
double dta_adt_compute_min(const DTA_SwitchedSystem* sys) {
    if (!sys || sys->n_modes <= 0) return INFINITY;
    int n = sys->state_dim, m = sys->n_modes;
    double tau_a_star = 0.0;
    int i;
    for (i = 0; i < m; i++) {
        double* Q = safe_alloc((size_t)(n*n) * sizeof(double));
        memset(Q, 0, (size_t)(n*n) * sizeof(double));
        int j;
        for (j = 0; j < n; j++) Q[j*n + j] = 1.0;
        double* P = safe_alloc((size_t)(n*n) * sizeof(double));
        int ret = dta_solve_lyapunov(sys->modes[i].A, n, Q, P);
        if (ret != 0) {
            free(Q); free(P);
            return INFINITY;
        }
        double* re = safe_alloc((size_t)n * sizeof(double));
        double* im = safe_alloc((size_t)n * sizeof(double));
        dta_eigenvalues(P, n, re, im, 200);
        double min_eig = re[0], max_eig = re[0];
        for (j = 1; j < n; j++) {
            if (re[j] < min_eig) min_eig = re[j];
            if (re[j] > max_eig) max_eig = re[j];
        }
        free(re); free(im);
        double cond = (min_eig > 1e-12) ? sqrt(max_eig / min_eig) : 1e6;
        double lambda = (min_eig > 1e-12) ? 1.0 / (2.0 * max_eig) : 0.0;
        double tau_i = (lambda > 0 && cond > 1.0) ? log(cond) / lambda : 0.0;
        if (tau_i > tau_a_star) tau_a_star = tau_i;
        free(Q); free(P);
    }
    return tau_a_star;
}

/* Mode-dependent minimum ADT */
void dta_adt_compute_mode_dependent(const DTA_SwitchedSystem* sys,
                                     double* tau_a_i_out) {
    if (!sys || !tau_a_i_out) return;
    int n = sys->state_dim, m = sys->n_modes;
    int i, j;
    for (i = 0; i < m; i++) tau_a_i_out[i] = 0.0;
    for (i = 0; i < m; i++) {
        double* Q = safe_alloc((size_t)(n*n) * sizeof(double));
        memset(Q, 0, (size_t)(n*n) * sizeof(double));
        for (j = 0; j < n; j++) Q[j*n + j] = 1.0;
        double* P = safe_alloc((size_t)(n*n) * sizeof(double));
        if (dta_solve_lyapunov(sys->modes[i].A, n, Q, P) == 0) {
            double* re = safe_alloc((size_t)n * sizeof(double));
            double* im = safe_alloc((size_t)n * sizeof(double));
            dta_eigenvalues(P, n, re, im, 200);
            double min_eig = re[0], max_eig = re[0];
            for (j = 1; j < n; j++) {
                if (re[j] < min_eig) min_eig = re[j];
                if (re[j] > max_eig) max_eig = re[j];
            }
            free(re); free(im);
            double cond = (min_eig > 1e-12) ? sqrt(max_eig / min_eig) : 1e6;
            double lambda = (min_eig > 1e-12) ? 1.0 / (2.0 * max_eig) : 0.0;
            tau_a_i_out[i] = (lambda > 0 && cond > 1.0) ? log(cond) / lambda : 0.0;
        }
        free(Q); free(P);
    }
}

/* Check if a given switching signal satisfies ADT condition */
bool dta_adt_signal_satisfies(const DTA_SwitchingSignal* sig,
                               double tau_a, double N0) {
    return dta_signal_validate_avg_dwell(sig, tau_a, N0);
}

/* ADT stability analysis */
DTA_ADTResult dta_adt_analyze(const DTA_SwitchedSystem* sys,
                               const DTA_SwitchingSignal* sig,
                               double N0) {
    DTA_ADTResult result;
    memset(&result, 0, sizeof(result));
    result.verdict = DTA_INCONCLUSIVE;
    result.N0 = N0;

    if (!sys || !sig) return result;

    result.actual_tau_a = dta_signal_avg_dwell_at(sig,
        sig->switch_times[sig->n_switches-1]);
    result.tau_a_min = dta_adt_compute_min(sys);

    DTA_DwellStabilityResult dwell_res = dta_analyze_dwell_stability(
        sys, result.actual_tau_a);
    result.verdict = dwell_res.verdict;
    result.decay_rate = dwell_res.decay_rate;
    if (dwell_res.lambda_i) free(dwell_res.lambda_i);
    if (dwell_res.common_P) free(dwell_res.common_P);

    result.mode_dependent_ok = (result.actual_tau_a > result.tau_a_min);
    return result;
}

/* Compute chatter bound N0 from a signal:
 * N0 = max_{T>=t>=0} { N_sigma(T,t) - (T-t)/tau_a } */
double dta_adt_compute_N0(const DTA_SwitchingSignal* sig, double tau_a) {
    if (!sig || sig->n_switches < 2 || tau_a <= 0) return 0.0;
    double N0 = 0.0;
    int i, j;
    for (i = 0; i < sig->n_switches; i++) {
        for (j = i + 1; j < sig->n_switches; j++) {
            double T = sig->switch_times[j];
            double t = sig->switch_times[i];
            if (T <= t) continue;
            int N = dta_signal_count_switches(sig, t, T);
            double bound = (double)N - (T - t) / tau_a;
            if (bound > N0) N0 = bound;
        }
    }
    return N0 > 0 ? N0 : 0.0;
}

/* Generate ADT-constrained switching signal */
DTA_SwitchingSignal* dta_adt_generate_signal(double t_start, double t_end,
    int n_modes, double tau_a, double N0) {
    if (n_modes <= 0 || tau_a <= 0 || t_end <= t_start) return NULL;
    return dta_signal_average_dwell(t_start, t_end, n_modes, tau_a, N0);
}

/* Convergence rate bound: alpha = lambda - ln(mu)/tau_a > 0 */
double dta_adt_convergence_rate(const DTA_SwitchedSystem* sys,
                                 double tau_a, double N0) {
    if (!sys || tau_a <= 0) return 0.0;
    DTA_DwellStabilityResult res = dta_analyze_dwell_stability(sys, tau_a);
    double rate = res.decay_rate;
    if (res.lambda_i) free(res.lambda_i);
    if (res.common_P) free(res.common_P);
    (void)N0;
    return rate;
}

/* Compare constant dwell vs ADT */
void dta_adt_compare(const DTA_SwitchingSignal* sig,
                      double* out_const_dwell, double* out_avg_dwell) {
    if (!sig) {
        if (out_const_dwell) *out_const_dwell = INFINITY;
        if (out_avg_dwell) *out_avg_dwell = INFINITY;
        return;
    }
    if (out_const_dwell)
        *out_const_dwell = dta_signal_min_dwell(sig);
    if (out_avg_dwell)
        *out_avg_dwell = dta_signal_avg_dwell_at(sig,
            sig->switch_times[sig->n_switches-1]);
}

void dta_adt_report(const DTA_ADTResult* result) {
    if (!result) return;
    printf("=== ADT Stability Report ===\n");
    printf("  Minimum ADT required: %.6f\n", result->tau_a_min);
    printf("  Actual ADT:           %.6f\n", result->actual_tau_a);
    printf("  Chatter bound N0:     %.6f\n", result->N0);
    printf("  Decay rate:           %.6f\n", result->decay_rate);
    printf("  Verdict:              %d\n", result->verdict);
    printf("  Mode-dependent OK:    %s\n",
           result->mode_dependent_ok ? "yes" : "no");
}
