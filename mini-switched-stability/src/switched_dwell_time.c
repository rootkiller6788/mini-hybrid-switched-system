#include "switched_types.h"
#include "switched_dwell_time.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* ============================================================================
 * L4: Minimum Dwell Time Computation
 * ============================================================================ */

double sdt_compute_min_dwell(double lambda_0, double mu) {
    if (lambda_0 <= 1e-12) return INFINITY;
    if (mu < 1.0) mu = 1.0;
    return log(mu) / lambda_0;
}

double sdt_compute_avg_dwell(double lambda_0, double mu) {
    return sdt_compute_min_dwell(lambda_0, mu);
}

double sdt_actual_avg_dwell(const SwitchingSignal *signal) {
    if (!signal || signal->n_switches <= 0) return INFINITY;
    double total = signal->switch_times[signal->n_switches] - signal->switch_times[0];
    if (total <= 1e-12) return 0.0;
    return total / (double)signal->n_switches;
}

int sdt_chatter_bound(const SwitchingSignal *signal) {
    if (!signal) return 0;
    /* N_0 is typically 0 or 1. For systems with maximum switching frequency,
     * compute the maximum number of switches in any interval of length tau_a */
    if (signal->n_switches <= 1) return 0;

    /* Estimate chatter: count max switches in any sliding window of size tau_avg */
    double tau_avg = sdt_actual_avg_dwell(signal);
    if (tau_avg <= 1e-12) return signal->n_switches;

    int max_chatter = 0;
    for (int k = 0; k < signal->n_switches; k++) {
        double t_start = signal->switch_times[k];
        double t_end = t_start + tau_avg;
        int count = 0;
        for (int j = k + 1; j <= signal->n_switches && signal->switch_times[j] <= t_end; j++) {
            count++;
        }
        if (count > max_chatter) max_chatter = count;
    }
    return max_chatter;
}

bool sdt_check_dwell(const double *switch_times, int n_switches, double tau_d) {
    if (!switch_times || n_switches <= 0) return true;
    for (int k = 0; k < n_switches; k++) {
        double interval = switch_times[k + 1] - switch_times[k];
        if (interval < tau_d - 1e-12) return false;
    }
    return true;
}

bool sdt_check_avg_dwell(const SwitchingSignal *signal, double tau_a, int N0) {
    if (!signal) return true;
    if (tau_a <= 1e-12) return false;

    /* Check N_sigma(T, t) <= N_0 + (T - t) / tau_a for all T > t >= 0 */
    for (int k = 0; k <= signal->n_switches; k++) {
        double T = signal->switch_times[k];
        int N_s = k; /* Number of switches in [0, T] */
        if (T > 0) {
            double bound = (double)N0 + T / tau_a;
            if ((double)N_s > bound + 1e-12) return false;
        }
    }
    return true;
}

/* ============================================================================
 * L4: Full Dwell-Time Analysis
 * ============================================================================ */

void sdt_full_analysis(SwitchedSystem *sys, DwellTimeAnalysis *dta) {
    if (!sys || !dta) return;

    /* Step 1: Compute stability margin lambda_0 */
    dta->stability_margin = sdt_stability_margin(sys);
    double lambda_0 = dta->stability_margin;

    /* Step 2: Compute Lyapunov functions and mismatch parameter mu */
    dta->mu_lyap = 1.0;

    if (sys->n_modes > 0) {
        MultipleLyapunovFunctions mlf;
        mlf.P = (SwitchedMatrix *)malloc((size_t)sys->n_modes * sizeof(SwitchedMatrix));
        mlf.min_eig = (double *)calloc((size_t)sys->n_modes, sizeof(double));
        mlf.max_eig = (double *)calloc((size_t)sys->n_modes, sizeof(double));
        mlf.is_valid = (bool *)calloc((size_t)sys->n_modes, sizeof(bool));
        mlf.decay_rates = (double *)calloc((size_t)sys->n_modes, sizeof(double));
        mlf.n_modes = sys->n_modes;
        mlf.mu = 1.0;
        mlf.mlf_condition = false;

        for (int i = 0; i < sys->n_modes; i++) {
            mlf.P[i] = sm_create(sys->state_dim, sys->state_dim);
        }

        /* Compute Lyapunov functions by solving A_i^T P_i + P_i A_i = -I */
        SwitchedMatrix Q = sm_create(sys->state_dim, sys->state_dim);
        sm_identity(&Q);
        for (int i = 0; i < sys->n_modes; i++) {
            if (sys->modes[i] && sys->modes[i]->is_hurwitz) {
                /* Need sss_solve_lyapunov from switched_stability.h, use inline version */
                bool solved = false;
                int n = sys->state_dim;
                if (n <= 4) {
                    /* Use Kronecker product method */
                    int n2 = n * n;
                    double *Kmat = (double *)calloc((size_t)(n2 * n2), sizeof(double));
                    double *rhs = (double *)malloc((size_t)n2 * sizeof(double));
                    if (Kmat && rhs) {
                        for (int ri = 0; ri < n; ri++) {
                            for (int rj = 0; rj < n; rj++) {
                                int row = ri * n + rj;
                                for (int ck = 0; ck < n; ck++)
                                    Kmat[row * n2 + ri * n + ck] += sys->modes[i]->A.data[ck * n + rj];
                                for (int ck = 0; ck < n; ck++)
                                    Kmat[row * n2 + ck * n + rj] += sys->modes[i]->A.data[ck * n + ri];
                                rhs[row] = -Q.data[ri * n + rj];
                            }
                        }
                        /* Simple Gauss-Seidel for solution */
                        double *x = (double *)calloc((size_t)n2, sizeof(double));
                        for (int gs_iter = 0; gs_iter < 1000; gs_iter++) {
                            double max_diff = 0.0;
                            for (int ri = 0; ri < n2; ri++) {
                                double old = x[ri];
                                double sum = rhs[ri];
                                for (int cj = 0; cj < n2; cj++)
                                    if (cj != ri) sum -= Kmat[ri * n2 + cj] * x[cj];
                                double diag = Kmat[ri * n2 + ri];
                                if (fabs(diag) > 1e-12) x[ri] = sum / diag;
                                double diff = fabs(x[ri] - old);
                                if (diff > max_diff) max_diff = diff;
                            }
                            if (max_diff < 1e-10) { solved = true; break; }
                        }
                        if (solved) {
                            for (int ri = 0; ri < n; ri++)
                                for (int rj = 0; rj < n; rj++)
                                    mlf.P[i].data[ri * n + rj] = x[ri * n + rj];
                            mlf.is_valid[i] = true;
                        }
                        free(x);
                        free(Kmat);
                        free(rhs);
                    }
                }
                if (solved) {
                    /* Estimate eigenvalue bounds */
                    double min_e = INFINITY, max_e = 0.0;
                    for (int r = 0; r < sys->state_dim; r++) {
                        double diag = mlf.P[i].data[r * sys->state_dim + r];
                        double radius = 0.0;
                        for (int c = 0; c < sys->state_dim; c++)
                            if (c != r) radius += fabs(mlf.P[i].data[r * sys->state_dim + c]);
                        if (diag - radius < min_e) min_e = diag - radius;
                        if (diag + radius > max_e) max_e = diag + radius;
                    }
                    mlf.min_eig[i] = (min_e > 0) ? min_e : 0.001;
                    mlf.max_eig[i] = max_e;
                    mlf.decay_rates[i] = 0.5 / mlf.max_eig[i];
                }
            }
        }
        sm_free(&Q);

        dta->mu_lyap = sdt_lyap_mismatch(&mlf);

        for (int i = 0; i < sys->n_modes; i++) sm_free(&mlf.P[i]);
        free(mlf.P); free(mlf.min_eig); free(mlf.max_eig);
        free(mlf.is_valid); free(mlf.decay_rates);
    }

    /* Step 3: Compute dwell time bounds */
    dta->tau_d = sdt_compute_min_dwell(lambda_0, dta->mu_lyap);
    dta->required_tau_a = sdt_compute_avg_dwell(lambda_0, dta->mu_lyap);
    dta->N0 = 1;

    /* Step 4: Evaluate actual switching signal */
    if (sys->signal && sys->signal->n_switches > 0) {
        dta->tau_a = sdt_actual_avg_dwell(sys->signal);
        dta->slow_enough = sdt_check_avg_dwell(sys->signal, dta->required_tau_a, dta->N0);
    }
}

double sdt_stability_margin(const SwitchedSystem *sys) {
    if (!sys || sys->n_modes <= 0) return 0.0;

    double lambda_0 = INFINITY;
    for (int i = 0; i < sys->n_modes; i++) {
        if (!sys->modes[i] || !sys->modes[i]->is_hurwitz) continue;
        int n = sys->state_dim;
        QRWorkspace *qr = qr_create(n);
        EigenvalueResult *eig = (EigenvalueResult *)malloc((size_t)n * sizeof(EigenvalueResult));
        qr_eigenvalues(qr, &sys->modes[i]->A, eig);

        double alpha_i = INFINITY;
        for (int j = 0; j < n; j++) {
            if (eig[j].real >= -1e-12) {
                alpha_i = -INFINITY;
                break;
            }
            double decay = -eig[j].real;
            if (decay < alpha_i) alpha_i = decay;
        }
        free(eig);
        qr_free(qr);

        if (alpha_i > 0 && alpha_i < lambda_0) lambda_0 = alpha_i;
    }

    return (lambda_0 > 0 && lambda_0 < INFINITY) ? lambda_0 : 0.0;
}

double sdt_lyap_mismatch(const MultipleLyapunovFunctions *mlf) {
    if (!mlf || mlf->n_modes <= 0) return 1.0;

    double mu = 1.0;
    for (int i = 0; i < mlf->n_modes; i++) {
        if (!mlf->is_valid[i]) continue;
        for (int j = 0; j < mlf->n_modes; j++) {
            if (!mlf->is_valid[j]) continue;
            if (mlf->min_eig[i] > 1e-12) {
                double ratio = mlf->max_eig[j] / mlf->min_eig[i];
                if (ratio > mu) mu = ratio;
            }
        }
    }
    return (mu < 1.0) ? 1.0 : mu;
}

SwitchedStabilityType sdt_certify(const DwellTimeAnalysis *dta, const SwitchingSignal *signal) {
    if (!dta || !signal) return SSTAB_UNSTABLE;

    /* Check minimum dwell time */
    bool min_dwell_ok = sdt_check_dwell(signal->switch_times, signal->n_switches, dta->tau_d);

    /* Check average dwell time */
    bool avg_dwell_ok = sdt_check_avg_dwell(signal, dta->required_tau_a, dta->N0);

    if (min_dwell_ok) {
        return SSTAB_GUES;
    } else if (avg_dwell_ok) {
        return SSTAB_EXPONENTIAL_STABLE;
    } else if (dta->stability_margin > 1e-12) {
        return SSTAB_ASYMPTOTIC_STABLE;
    }

    return SSTAB_UNSTABLE;
}

double sdt_critical_dwell_time(SwitchedSystem *sys, double T_max) {
    if (!sys || sys->n_modes <= 0) return INFINITY;

    DwellTimeAnalysis dta;
    sdt_full_analysis(sys, &dta);

    if (dta.stability_margin <= 1e-12) return INFINITY;

    /* Binary search for critical dwell time */
    double lo = 0.0, hi = T_max;
    for (int iter = 0; iter < 50; iter++) {
        double mid = (lo + hi) / 2.0;
        /* Check if avg dwell time condition is sufficient */
        if (mid > dta.required_tau_a) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    return hi;
}

double sdt_dwell_for_tolerance(double lambda_0, double eps, double mu) {
    if (lambda_0 <= 1e-12 || eps <= 1e-12) return INFINITY;
    /* ||x(T)|| <= exp(-lambda_0 * T + N_0 * ln(mu)) * ||x(0)|| */
    /* Want: exp(-lambda_0 * T + ln(mu) * T/tau_d) <= eps */
    /* Solve for tau_d: */
    /* -lambda_0 * tau_d + ln(mu) <= ln(eps)/N_s */
    /* Conservative: tau_d >= ln(mu) / lambda_0 + |ln(eps)|/lambda_0 */
    double base_tau = log(mu) / lambda_0;
    double extra = -log(eps) / lambda_0;
    return base_tau + extra;
}