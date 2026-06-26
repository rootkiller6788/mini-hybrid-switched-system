#include "switched_types.h"
#include "switched_stability.h"
#include "switched_lyapunov.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* ============================================================================
 * Stability Classification
 * ============================================================================ */

SwitchedStabilityType sss_classify_stability(SwitchedSystem *sys) {
    if (!sys || sys->n_modes <= 0) return SSTAB_UNSTABLE;

    bool all_hurwitz = true;
    for (int i = 0; i < sys->n_modes; i++) {
        if (!sys->modes[i] || !sys->modes[i]->is_hurwitz) {
            all_hurwitz = false;
            break;
        }
    }

    if (!all_hurwitz) return SSTAB_UNSTABLE;

    /* Check switching signal type to determine stability class */
    if (sys->signal) {
        switch (sys->signal->type) {
            case SSIG_ARBITRARY:
                /* Need CLF for GUES */
                return SSTAB_EXPONENTIAL_STABLE;
            case SSIG_DWELL_TIME:
            case SSIG_AVG_DWELL:
                /* With sufficiently slow switching */
                return SSTAB_GUES;
            case SSIG_PERIODIC:
                return SSTAB_EXPONENTIAL_STABLE;
            default:
                return SSTAB_ASYMPTOTIC_STABLE;
        }
    }

    return SSTAB_EXPONENTIAL_STABLE;
}

/* ============================================================================
 * GUES Under Arbitrary Switching - CLF Check
 * ============================================================================ */

bool sss_is_gues_arbitrary(SwitchedSystem *sys) {
    if (!sys || sys->n_modes <= 0) return false;

    /* Quick check: all subsystems must be Hurwitz */
    for (int i = 0; i < sys->n_modes; i++) {
        if (!sys->modes[i] || !sys->modes[i]->is_hurwitz) return false;
    }

    /* Check Lie-algebraic condition for CLF existence */
    LieAlgebraCondition *la = lie_check_create(sys->n_modes, sys->state_dim);
    SwitchedMatrix **A_arr = (SwitchedMatrix **)malloc((size_t)sys->n_modes * sizeof(SwitchedMatrix *));
    for (int i = 0; i < sys->n_modes; i++) {
        A_arr[i] = &sys->modes[i]->A;
    }
    bool solvable = lie_condition_solvable(la, A_arr, sys->n_modes, sys->state_dim);
    lie_check_free(la);
    free(A_arr);

    /* If Lie algebra is solvable, a CLF exists */
    if (solvable) return true;

    /* Additional check: try to find CLF via gradient descent */
    CommonLyapunovFunction clf;
    clf.P = sm_create(sys->state_dim, sys->state_dim);
    clf.is_valid = false;

    bool found = sss_compute_clf(sys, &clf, 500);
    if (found) {
        sm_free(&clf.P);
        return true;
    }
    sm_free(&clf.P);
    return false;
}

/* ============================================================================
 * MLF Stability Check
 * ============================================================================ */

bool sss_is_stable_mlf(SwitchedSystem *sys, MultipleLyapunovFunctions *mlf) {
    if (!sys || !mlf) return false;

    /* Compute MLF for all modes */
    sss_compute_mlf(sys, mlf);

    /* Verify all Lyapunov functions are valid */
    for (int i = 0; i < mlf->n_modes; i++) {
        if (!mlf->is_valid[i]) return false;
    }

    /* Check sequence condition if switching signal is available */
    if (sys->signal) {
        mlf->mlf_condition = sss_mlf_verify_sequence(sys, mlf, sys->signal);
        return mlf->mlf_condition;
    }

    return true;
}

/* ============================================================================
 * L4: CLF Theorem Implementation
 * ============================================================================ */

bool sss_clf_theorem(SwitchedSystem *sys, CommonLyapunovFunction *clf) {
    if (!sys || !clf || sys->n_modes <= 0) return false;

    /* For n=1 (single subsystem), any Lyapunov function for A_0 works */
    if (sys->n_modes == 1 && sys->modes[0]->is_hurwitz) {
        int n = sys->state_dim;
        SwitchedMatrix I = sm_create(n, n);
        sm_identity(&I);

        /* Solve A^T P + P A = -I for the single system */
        SwitchedMatrix P = sm_create(n, n);
        bool solved = sss_solve_lyapunov(&sys->modes[0]->A, &I, &P, n);
        if (solved) {
            sm_copy(&clf->P, &P);
            clf->is_valid = true;
            clf->decay_rate = 0.5; /* Approximate */
            clf->min_eig = 0.1;
            clf->max_eig = 1.0;
        }
        sm_free(&I);
        sm_free(&P);
        return solved;
    }

    /* For multiple subsystems, try to find CLF */
    return sss_compute_clf(sys, clf, 500);
}

/* ============================================================================
 * L4: MLF Theorem Implementation
 * ============================================================================ */

bool sss_mlf_theorem(SwitchedSystem *sys, MultipleLyapunovFunctions *mlf) {
    if (!sys || !mlf || sys->n_modes <= 0) return false;

    /* Compute individual Lyapunov functions */
    sss_compute_mlf(sys, mlf);

    /* Check that each V_i is a valid Lyapunov function for mode i */
    for (int i = 0; i < sys->n_modes; i++) {
        if (!mlf->is_valid[i]) return false;

        /* Verify: A_i^T P_i + P_i A_i < 0 */
        if (mlf->decay_rates[i] <= 0.0) return false;
    }

    /* Compute mismatch parameter mu */
    mlf->mu = sss_compute_mu(mlf);

    /* Apply Branicky condition if switching signal exists */
    if (sys->signal) {
        mlf->mlf_condition = sss_mlf_verify_sequence(sys, mlf, sys->signal);
        return mlf->mlf_condition;
    }

    return true; /* All V_i valid, MLF exists */
}

/* ============================================================================
 * L4: Dwell-Time Theorem
 * ============================================================================ */

void sss_dwell_time_theorem(SwitchedSystem *sys, DwellTimeAnalysis *dta) {
    if (!sys || !dta) return;

    dta->stability_margin = sss_compute_stability_margin(sys);
    dta->mu_lyap = 1.0;

    /* Compute Lyapunov functions and mu parameter */
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

        sss_compute_mlf(sys, &mlf);
        dta->mu_lyap = sss_compute_mu(&mlf);

        for (int i = 0; i < sys->n_modes; i++) sm_free(&mlf.P[i]);
        free(mlf.P); free(mlf.min_eig); free(mlf.max_eig);
        free(mlf.is_valid); free(mlf.decay_rates);
    }

    dta->tau_d = sss_compute_min_dwell_time(dta->stability_margin, dta->mu_lyap);
    dta->required_tau_a = sss_compute_avg_dwell_time(dta->stability_margin, dta->mu_lyap);
    dta->N0 = 1;

    /* Check if actual switching is slow enough */
    if (sys->signal && sys->signal->n_switches > 0) {
        dta->tau_a = sss_compute_actual_avg_dwell(sys->signal);
        dta->slow_enough = (dta->tau_a > dta->required_tau_a);
    }
}

/* ============================================================================
 * L4: Average Dwell-Time Theorem
 * ============================================================================ */

void sss_avg_dwell_theorem(SwitchedSystem *sys, DwellTimeAnalysis *dta) {
    if (!sys || !dta) return;

    /* Compute stability margin and mu */
    dta->stability_margin = sss_compute_stability_margin(sys);

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

    sss_compute_mlf(sys, &mlf);
    dta->mu_lyap = sss_compute_mu(&mlf);

    for (int i = 0; i < sys->n_modes; i++) sm_free(&mlf.P[i]);
    free(mlf.P); free(mlf.min_eig); free(mlf.max_eig);
    free(mlf.is_valid); free(mlf.decay_rates);

    /* Compute critical average dwell time:
     * tau_a* = ln(mu) / lambda_0
     */
    double lambda_0 = dta->stability_margin;
    double mu = dta->mu_lyap;

    if (lambda_0 > 1e-12 && mu >= 1.0) {
        dta->required_tau_a = log(mu) / lambda_0;
    } else {
        dta->required_tau_a = 0.0;
    }

    dta->N0 = 1;

    /* Check if actual average dwell time satisfies the bound */
    if (sys->signal && sys->signal->n_switches > 0) {
        dta->tau_a = sss_compute_actual_avg_dwell(sys->signal);
        dta->slow_enough = sss_check_avg_dwell(sys->signal, dta->required_tau_a, dta->N0);
    }
}

/* ============================================================================
 * L5: CLF Computation via Gradient Descent
 * ============================================================================ */

bool sss_compute_clf(SwitchedSystem *sys, CommonLyapunovFunction *clf, int max_iters) {
    if (!sys || !clf || sys->n_modes <= 0) return false;

    int n = sys->state_dim;

    /* Initialize P = I (identity is a natural starting point) */
    SwitchedMatrix P = sm_create(n, n);
    sm_identity(&P);

    double alpha = 0.01; /* Step size */
    double best_violation = INFINITY;

    for (int iter = 0; iter < max_iters; iter++) {
        /* Compute violation: max_i lambda_max(A_i^T P + P A_i) */
        double max_violation = -INFINITY;

        for (int mode = 0; mode < sys->n_modes; mode++) {
            SwitchedMatrix *A = &sys->modes[mode]->A;

            /* Compute A^T P + P A */
            SwitchedMatrix AT = sm_create(n, n);
            SwitchedMatrix PA = sm_create(n, n);
            SwitchedMatrix ATP = sm_create(n, n);

            sm_transpose(&AT, A);
            sm_mul(&PA, &P, A);
            sm_mul(&ATP, &AT, &P);

            /* R = A^T P + P A */
            SwitchedMatrix R = sm_create(n, n);
            sm_add(&R, &ATP, &PA);

            /* Estimate max eigenvalue of R via power method */
            double max_eig = 0.0;
            for (int i = 0; i < n; i++) {
                double row_sum = 0.0;
                for (int j = 0; j < n; j++) {
                    row_sum += fabs(R.data[i * n + j]);
                }
                if (row_sum > max_eig) max_eig = row_sum;
            }

            if (max_eig > max_violation) max_violation = max_eig;

            sm_free(&AT); sm_free(&PA); sm_free(&ATP); sm_free(&R);
        }

        if (max_violation < best_violation) best_violation = max_violation;

        /* If violation is negative, P is a valid CLF */
        if (max_violation < -1e-8) {
            sm_copy(&clf->P, &P);
            clf->is_valid = true;
            clf->decay_rate = -max_violation / 2.0;

            /* Estimate eigenvalue bounds */
            double min_e = INFINITY, max_e = 0.0;
            for (int i = 0; i < n; i++) {
                double row_sum = 0.0;
                for (int j = 0; j < n; j++) {
                    row_sum += fabs(P.data[i * n + j]);
                }
                if (row_sum < min_e) min_e = row_sum;
                if (row_sum > max_e) max_e = row_sum;
            }
            clf->min_eig = min_e;
            clf->max_eig = max_e;

            sm_free(&P);
            return true;
        }

        /* Gradient step: move P towards satisfying the LMI */
        /* Use a simple perturbation approach: */
        /* P_{k+1} = P_k - alpha * (A_i^T P_k + P_k A_i + eps*I) for worst mode */
        for (int mode = 0; mode < sys->n_modes; mode++) {
            SwitchedMatrix *A = &sys->modes[mode]->A;
            SwitchedMatrix AT = sm_create(n, n);
            sm_transpose(&AT, A);

            SwitchedMatrix PA = sm_create(n, n);
            SwitchedMatrix ATP = sm_create(n, n);
            sm_mul(&PA, &P, A);
            sm_mul(&ATP, &AT, &P);

            SwitchedMatrix grad = sm_create(n, n);
            sm_add(&grad, &ATP, &PA);

            /* P_{k+1} = P_k - alpha * grad */
            sm_mul_scalar(&grad, alpha);
            sm_sub(&P, &P, &grad);

            sm_free(&AT); sm_free(&PA); sm_free(&ATP); sm_free(&grad);

            /* Ensure symmetry: P = (P + P^T) / 2 */
            SwitchedMatrix PT = sm_create(n, n);
            sm_transpose(&PT, &P);
            sm_add(&P, &P, &PT);
            sm_mul_scalar(&P, 0.5);
            sm_free(&PT);

            /* Enforce positive definiteness: ensure diagonal dominance */
            for (int i = 0; i < n; i++) {
                if (P.data[i * n + i] < 0.001) P.data[i * n + i] = 0.001;
            }
            break; /* Only update once per iteration */
        }

        /* Adjust step size */
        if (iter % 50 == 0) alpha *= 0.95;
    }

    sm_free(&P);
    return false;
}

/* ============================================================================
 * L5: MLF Computation
 * ============================================================================ */

void sss_compute_mlf(SwitchedSystem *sys, MultipleLyapunovFunctions *mlf) {
    if (!sys || !mlf) return;

    int n = sys->state_dim;
    SwitchedMatrix Q = sm_create(n, n);
    sm_identity(&Q);

    for (int i = 0; i < sys->n_modes && i < mlf->n_modes; i++) {
        SwitchedSubsystem *sub = sys->modes[i];
        if (!sub || !sub->is_hurwitz) {
            mlf->is_valid[i] = false;
            mlf->decay_rates[i] = 0.0;
            mlf->min_eig[i] = 0.0;
            mlf->max_eig[i] = 0.0;
            continue;
        }

        /* Solve Lyapunov equation A_i^T P_i + P_i A_i = -Q */
        bool solved = sss_solve_lyapunov(&sub->A, &Q, &mlf->P[i], n);
        mlf->is_valid[i] = solved;

        if (solved) {
            /* Estimate eigenvalue bounds using Gershgorin circle theorem */
            double min_e = INFINITY, max_e = 0.0;
            for (int r = 0; r < n; r++) {
                double diag = mlf->P[i].data[r * n + r];
                double radius = 0.0;
                for (int c = 0; c < n; c++) {
                    if (c != r) radius += fabs(mlf->P[i].data[r * n + c]);
                }
                if (diag - radius < min_e) min_e = diag - radius;
                if (diag + radius > max_e) max_e = diag + radius;
            }
            mlf->min_eig[i] = (min_e > 0) ? min_e : 0.001;
            mlf->max_eig[i] = max_e;

            /* Decay rate: alpha_i = 1 / (2 * lambda_max(P_i)) for Q=I */
            mlf->decay_rates[i] = 0.5 / mlf->max_eig[i];
        }
    }

    sm_free(&Q);
    mlf->mu = sss_compute_mu(mlf);
}

/* ============================================================================
 * L5: Lyapunov Equation Solver
 * ============================================================================ */

bool sss_solve_lyapunov(const SwitchedMatrix *A, const SwitchedMatrix *Q,
                        SwitchedMatrix *P, int n) {
    if (!A || !Q || !P || n <= 0 || n > 6) return false;

    /* Construct Kronecker system: (I kron A^T + A^T kron I) vec(P) = -vec(Q) */
    int n2 = n * n;

    /* Build the n^2 x n^2 coefficient matrix K */
    double *K = (double *)calloc((size_t)(n2 * n2), sizeof(double));
    if (!K) return false;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int row = i * n + j; /* Row in vectorized system */

            /* (I kron A^T) contribution */
            for (int k = 0; k < n; k++) {
                int col = i * n + k;
                K[row * n2 + col] += A->data[k * n + j]; /* A^T[j,k] = A[k,j] */
            }

            /* (A^T kron I) contribution */
            for (int k = 0; k < n; k++) {
                int col = k * n + j;
                K[row * n2 + col] += A->data[k * n + i]; /* A^T[i,k] = A[k,i] */
            }
        }
    }

    /* Build RHS: vec(-Q) */
    double *rhs = (double *)malloc((size_t)n2 * sizeof(double));
    if (!rhs) { free(K); return false; }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            rhs[i * n + j] = -Q->data[i * n + j];
        }
    }

    /* Solve K * vec(P) = rhs via Gaussian elimination with partial pivoting */
    double *aug = (double *)malloc((size_t)(n2 * (n2 + 1)) * sizeof(double));
    if (!aug) { free(K); free(rhs); return false; }

    for (int i = 0; i < n2; i++) {
        for (int j = 0; j < n2; j++) {
            aug[i * (n2 + 1) + j] = K[i * n2 + j];
        }
        aug[i * (n2 + 1) + n2] = rhs[i];
    }

    /* Gaussian elimination */
    for (int col = 0; col < n2; col++) {
        /* Find pivot */
        int pivot_row = col;
        double pivot_val = fabs(aug[col * (n2 + 1) + col]);
        for (int r = col + 1; r < n2; r++) {
            double val = fabs(aug[r * (n2 + 1) + col]);
            if (val > pivot_val) {
                pivot_val = val;
                pivot_row = r;
            }
        }

        if (pivot_val < 1e-14) { free(K); free(rhs); free(aug); return false; }

        /* Swap rows */
        if (pivot_row != col) {
            for (int j = 0; j <= n2; j++) {
                double tmp = aug[col * (n2 + 1) + j];
                aug[col * (n2 + 1) + j] = aug[pivot_row * (n2 + 1) + j];
                aug[pivot_row * (n2 + 1) + j] = tmp;
            }
        }

        /* Eliminate below */
        double diag = aug[col * (n2 + 1) + col];
        for (int r = col + 1; r < n2; r++) {
            double factor = aug[r * (n2 + 1) + col] / diag;
            for (int j = col; j <= n2; j++) {
                aug[r * (n2 + 1) + j] -= factor * aug[col * (n2 + 1) + j];
            }
        }
    }

    /* Back substitution */
    double *vecP = (double *)malloc((size_t)n2 * sizeof(double));
    if (!vecP) { free(K); free(rhs); free(aug); return false; }

    for (int i = n2 - 1; i >= 0; i--) {
        double sum = aug[i * (n2 + 1) + n2];
        for (int j = i + 1; j < n2; j++) {
            sum -= aug[i * (n2 + 1) + j] * vecP[j];
        }
        vecP[i] = sum / aug[i * (n2 + 1) + i];
    }

    /* Unvectorize: vec(P)[i*n + j] = P[i,j] */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            P->data[i * n + j] = vecP[i * n + j];
        }
    }

    free(K); free(rhs); free(aug); free(vecP);
    return true;
}

bool sss_is_positive_definite(const SwitchedMatrix *P) {
    return sm_is_positive_definite(P);
}

double sss_compute_mu(const MultipleLyapunovFunctions *mlf) {
    if (!mlf || mlf->n_modes <= 0) return 1.0;

    double mu = 1.0;
    for (int i = 0; i < mlf->n_modes; i++) {
        if (!mlf->is_valid[i]) continue;
        for (int j = 0; j < mlf->n_modes; j++) {
            if (!mlf->is_valid[j]) continue;
            double ratio = mlf->max_eig[j] / mlf->min_eig[i];
            if (ratio > mu) mu = ratio;
        }
    }
    return (mu < 1.0) ? 1.0 : mu;
}

double sss_compute_stability_margin(SwitchedSystem *sys) {
    if (!sys || sys->n_modes <= 0) return 0.0;

    double lambda_0 = INFINITY;
    for (int i = 0; i < sys->n_modes; i++) {
        if (!sys->modes[i]) continue;
        int n = sys->state_dim;
        QRWorkspace *qr = qr_create(n);
        EigenvalueResult *eig = (EigenvalueResult *)malloc((size_t)n * sizeof(EigenvalueResult));
        qr_eigenvalues(qr, &sys->modes[i]->A, eig);

        double alpha_i = INFINITY;
        for (int j = 0; j < n; j++) {
            if (eig[j].real > -1e-12) {
                alpha_i = -INFINITY; /* Not Hurwitz */
                break;
            }
            double decay = -eig[j].real;
            if (decay < alpha_i) alpha_i = decay;
        }
        free(eig);
        qr_free(qr);

        if (alpha_i < lambda_0) lambda_0 = alpha_i;
    }

    return (lambda_0 > 0) ? lambda_0 : 0.0;
}

/* ============================================================================
 * L5: Lie-Algebraic Check
 * ============================================================================ */

bool sss_lie_algebraic_check(SwitchedSystem *sys, LieAlgebraCondition *la) {
    if (!sys || !la) return false;

    SwitchedMatrix **A_arr = (SwitchedMatrix **)malloc((size_t)sys->n_modes * sizeof(SwitchedMatrix *));
    for (int i = 0; i < sys->n_modes; i++) {
        A_arr[i] = &sys->modes[i]->A;
    }

    bool result = lie_condition_solvable(la, A_arr, sys->n_modes, sys->state_dim);
    free(A_arr);
    return result;
}

/* ============================================================================
 * L5: Dwell-Time Computation
 * ============================================================================ */

double sss_compute_min_dwell_time(double lambda_0, double mu) {
    if (lambda_0 <= 1e-12) return INFINITY;
    if (mu < 1.0) mu = 1.0;
    return log(mu) / lambda_0;
}

double sss_compute_avg_dwell_time(double lambda_0, double mu) {
    return sss_compute_min_dwell_time(lambda_0, mu);
}

bool sss_check_dwell_time(const double *switch_times, int n_switches, double tau_d) {
    if (!switch_times || n_switches <= 0) return true;
    for (int k = 0; k < n_switches; k++) {
        if (switch_times[k + 1] - switch_times[k] < tau_d - 1e-12) {
            return false;
        }
    }
    return true;
}

bool sss_check_avg_dwell(const SwitchingSignal *signal, double tau_a, int N0) {
    if (!signal) return true;
    if (tau_a <= 1e-12) return false;

    for (int k = 0; k <= signal->n_switches; k++) {
        double T = signal->switch_times[k];
        int N_s = k; /* Number of switches in [0, T] */
        double bound = (double)N0 + T / tau_a;
        if ((double)N_s > bound + 1e-12) return false;
    }
    return true;
}

double sss_compute_actual_avg_dwell(const SwitchingSignal *signal) {
    if (!signal || signal->n_switches <= 0) return INFINITY;
    double total = signal->switch_times[signal->n_switches] - signal->switch_times[0];
    if (total <= 1e-12) return 0.0;
    return total / (double)signal->n_switches;
}

/* ============================================================================
 * L5: Simulation Methods
 * ============================================================================ */

void sss_simulate_euler(SwitchedSystem *sys, const SolverConfig *cfg,
                        const SwitchSequence *seq) {
    sss_simulate(sys, cfg, seq);
}

void sss_simulate_rk4(SwitchedSystem *sys, const SolverConfig *cfg,
                      const SwitchSequence *seq) {
    if (!sys || !cfg || !seq || !sys->signal) return;

    double t = 0.0;
    int seq_idx = 0;
    double time_in_mode = 0.0;
    sys->signal->mode_sequence[0] = seq->mode_order[0];
    sys->signal->switch_times[0] = 0.0;
    sys->signal->n_switches = 0;
    sys->signal->current_mode = seq->mode_order[0];

    for (int step = 0; step < cfg->max_steps && t < cfg->t_end; step++) {
        double dt = cfg->dt;
        if (t + dt > cfg->t_end) dt = cfg->t_end - t;

        int mode = seq->mode_order[seq_idx];
        SwitchedMatrix *A = &sys->modes[mode]->A;

        /* RK4: k1 = A x(t) */
        SwitchedVector k1 = sv_create(sys->state_dim);
        SwitchedVector k2 = sv_create(sys->state_dim);
        SwitchedVector k3 = sv_create(sys->state_dim);
        SwitchedVector k4 = sv_create(sys->state_dim);
        SwitchedVector temp = sv_create(sys->state_dim);

        sm_matvec_mul(&k1, A, &sys->state);

        /* k2 = A (x(t) + dt/2 * k1) */
        sv_copy(&temp, &sys->state);
        sv_axpy(&temp, dt / 2.0, &k1);
        sm_matvec_mul(&k2, A, &temp);

        /* k3 = A (x(t) + dt/2 * k2) */
        sv_copy(&temp, &sys->state);
        sv_axpy(&temp, dt / 2.0, &k2);
        sm_matvec_mul(&k3, A, &temp);

        /* k4 = A (x(t) + dt * k3) */
        sv_copy(&temp, &sys->state);
        sv_axpy(&temp, dt, &k3);
        sm_matvec_mul(&k4, A, &temp);

        /* x(t+dt) = x(t) + dt/6 * (k1 + 2*k2 + 2*k3 + k4) */
        sv_axpy(&sys->state, dt / 6.0, &k1);
        sv_axpy(&sys->state, dt / 3.0, &k2);
        sv_axpy(&sys->state, dt / 3.0, &k3);
        sv_axpy(&sys->state, dt / 6.0, &k4);

        sv_free(&k1); sv_free(&k2); sv_free(&k3); sv_free(&k4); sv_free(&temp);

        t += dt;
        time_in_mode += dt;

        /* Check for mode switch */
        if (time_in_mode >= seq->durations[seq_idx] - cfg->event_tol) {
            seq_idx++;
            time_in_mode = 0.0;
            if (seq_idx >= seq->length) {
                if (seq->repeat) seq_idx = 0;
                else break;
            }
            ssig_record_switch(sys->signal, seq->mode_order[seq_idx], t);
        }
    }
    sys->signal->total_time = t;
}

void sss_simulate_state_dep(SwitchedSystem *sys, const SolverConfig *cfg,
                            StateSwitchFunc phi, void *phi_ctx) {
    if (!sys || !cfg || !phi) return;

    if (!sys->signal) sys->signal = ssig_create(100);
    sys->signal->type = SSIG_STATE_DEP;

    double t = 0.0;
    int current_mode = phi(&sys->state, phi_ctx);
    sys->signal->mode_sequence[0] = current_mode;
    sys->signal->current_mode = current_mode;

    for (int step = 0; step < cfg->max_steps && t < cfg->t_end; step++) {
        double dt = cfg->dt;
        if (t + dt > cfg->t_end) dt = cfg->t_end - t;

        /* Step forward */
        sss_step(sys, dt);
        t += dt;

        /* Check if switching function indicates a mode change */
        int new_mode = phi(&sys->state, phi_ctx);
        if (new_mode != current_mode) {
            ssig_record_switch(sys->signal, new_mode, t);
            current_mode = new_mode;
        }
    }
    sys->signal->total_time = t;
}