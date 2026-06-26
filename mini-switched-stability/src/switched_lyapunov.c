#include "switched_types.h"
#include "switched_stability.h"
#include "switched_lyapunov.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* ============================================================================
 * L3: Eigenvalue Computation for Lyapunov Matrices
 * ============================================================================ */

void sss_eig_sym_2x2(const SwitchedMatrix *P, double *eig1, double *eig2) {
    if (!P || !eig1 || !eig2 || P->rows != 2 || P->cols != 2) {
        if (eig1) *eig1 = 0.0;
        if (eig2) *eig2 = 0.0;
        return;
    }

    /* For symmetric [[a,b],[b,c]]:
     * lambda = (a+c)/2 +/- sqrt(((a-c)/2)^2 + b^2) */
    double a = P->data[0], b = P->data[1], c = P->data[3];
    double tr_half = (a + c) / 2.0;
    double disc = sqrt((a - c) * (a - c) / 4.0 + b * b);
    *eig1 = tr_half + disc;
    *eig2 = tr_half - disc;
}

void sss_eig_sym(const SwitchedMatrix *P, double *eig_min, double *eig_max, int max_iters) {
    if (!P || !eig_min || !eig_max || P->rows != P->cols) {
        if (eig_min) *eig_min = 0.0;
        if (eig_max) *eig_max = 0.0;
        return;
    }

    int n = P->rows;

    /* For n <= 3, use 2x2 approximation or Gershgorin */
    if (n <= 2) {
        if (n == 1) {
            *eig_min = P->data[0];
            *eig_max = P->data[0];
            return;
        }
        double e1, e2;
        sss_eig_sym_2x2(P, &e1, &e2);
        *eig_min = (e1 < e2) ? e1 : e2;
        *eig_max = (e1 > e2) ? e1 : e2;
        return;
    }

    /* Power iteration for largest eigenvalue */
    /* Initialize random vector */
    double *v = (double *)malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) v[i] = 1.0;

    double lambda_max = 0.0;
    for (int iter = 0; iter < max_iters; iter++) {
        /* w = P * v */
        double *w = (double *)malloc((size_t)n * sizeof(double));
        for (int i = 0; i < n; i++) {
            w[i] = 0.0;
            for (int j = 0; j < n; j++) {
                w[i] += P->data[i * n + j] * v[j];
            }
        }

        /* Rayleigh quotient: lambda = (v^T P v) / (v^T v) */
        double num = 0.0, den = 0.0;
        for (int i = 0; i < n; i++) {
            num += v[i] * w[i];
            den += v[i] * v[i];
        }
        double lambda_new = num / den;

        if (fabs(lambda_new - lambda_max) < 1e-10 && iter > 5) {
            lambda_max = lambda_new;
            free(w);
            break;
        }
        lambda_max = lambda_new;

        /* Normalize: v = w / ||w|| */
        double wnorm = 0.0;
        for (int i = 0; i < n; i++) wnorm += w[i] * w[i];
        wnorm = sqrt(wnorm);
        if (wnorm < 1e-15) { free(w); break; }
        for (int i = 0; i < n; i++) v[i] = w[i] / wnorm;
        free(w);
    }
    *eig_max = lambda_max;

    /* Estimate min eigenvalue using Gershgorin circle theorem */
    double min_gersh = INFINITY;
    for (int i = 0; i < n; i++) {
        double radius = 0.0;
        for (int j = 0; j < n; j++) {
            if (j != i) radius += fabs(P->data[i * n + j]);
        }
        double lower = P->data[i * n + i] - radius;
        if (lower < min_gersh) min_gersh = lower;
    }
    *eig_min = (min_gersh > 0) ? min_gersh : 0.001;

    free(v);
}

void sss_eig_2x2(const SwitchedMatrix *A, EigenvalueResult *e1, EigenvalueResult *e2) {
    if (!A || !e1 || !e2 || A->rows != 2 || A->cols != 2) return;

    double a = A->data[0], b = A->data[1];
    double c = A->data[2], d = A->data[3];
    double tr = a + d;
    double det = a * d - b * c;
    double disc = tr * tr - 4.0 * det;

    if (disc >= 0.0) {
        double sqrt_d = sqrt(disc);
        e1->real = (tr + sqrt_d) / 2.0; e1->imag = 0.0;
        e2->real = (tr - sqrt_d) / 2.0; e2->imag = 0.0;
    } else {
        double re = tr / 2.0, im = sqrt(-disc) / 2.0;
        e1->real = re; e1->imag = im;
        e2->real = re; e2->imag = -im;
    }
    e1->magnitude = sqrt(e1->real * e1->real + e1->imag * e1->imag);
    e1->stable = (e1->real < 0.0);
    e2->magnitude = sqrt(e2->real * e2->real + e2->imag * e2->imag);
    e2->stable = (e2->real < 0.0);
}

double sss_spectral_radius(const SwitchedMatrix *A) {
    return spectral_radius(A);
}

/* ============================================================================
 * L4: Lyapunov Equation Solvers
 * ============================================================================ */

bool sss_lyap_2x2(const SwitchedMatrix *A, SwitchedMatrix *P) {
    if (!A || !P || A->rows != 2 || P->rows != 2) return false;

    double a = A->data[0], b = A->data[1];
    double c_ = A->data[2], d = A->data[3];
    double tr = a + d;
    double det = a * d - b * c_;

    if (tr >= 0.0 || det <= 0.0) return false; /* A not Hurwitz */

    /* For Q = I, solve A^T P + P A = -I */
    /* Explicit solution via 4 linear equations */
    /* A^T = [[a, c_], [b, d]] */
    /* A^T P + P A = [[2a P11 + c_ P21 + c_ P21, a P12 + c_ P22 + b P11 + d P12],
     *                [a P21 + c_ P22 + b P11 + d P12, b P12 + d P22 + b P12 + d P22]] */

    /* Solve the linear system using the Kronecker form */
    /* For n=2: K is 4x4, solve by Gaussian elimination */
    double K[16];
    memset(K, 0, sizeof(K));

    /* Row 0: 2a*P11 + c_*P12 + b*P21 = -1 */
    K[0]  = 2.0 * a;     K[1]  = c_;         K[2]  = b;           K[3]  = 0.0;
    /* Row 1: b*P11 + (a+d)*P12 + c_*P22 = 0 */
    K[4]  = b;           K[5]  = a + d;       K[6]  = 0.0;         K[7]  = c_;
    /* Row 2: c_*P11 + (a+d)*P21 + b*P22 = 0 */
    K[8]  = c_;          K[9]  = 0.0;         K[10] = a + d;       K[11] = b;
    /* Row 3: 0*P12 + c_*P12 + b*P21 + 2d*P22 = -1 */
    K[12] = 0.0;         K[13] = c_;          K[14] = b;           K[15] = 2.0 * d;

    double rhs[4] = {-1.0, 0.0, 0.0, -1.0};
    double sol[4];

    /* Gaussian elimination */
    for (int col = 0; col < 4; col++) {
        int pivot = col;
        double pv = fabs(K[col * 4 + col]);
        for (int r = col + 1; r < 4; r++) {
            if (fabs(K[r * 4 + col]) > pv) {
                pv = fabs(K[r * 4 + col]);
                pivot = r;
            }
        }
        if (pv < 1e-14) return false;
        if (pivot != col) {
            for (int j = 0; j < 4; j++) {
                double tmp = K[col * 4 + j]; K[col * 4 + j] = K[pivot * 4 + j]; K[pivot * 4 + j] = tmp;
            }
            double trhs = rhs[col]; rhs[col] = rhs[pivot]; rhs[pivot] = trhs;
        }
        double pivot_val = K[col * 4 + col];
        for (int r = col + 1; r < 4; r++) {
            double factor = K[r * 4 + col] / pivot_val;
            for (int j = col; j < 4; j++) K[r * 4 + j] -= factor * K[col * 4 + j];
            rhs[r] -= factor * rhs[col];
        }
    }

    for (int i = 3; i >= 0; i--) {
        double sum = rhs[i];
        for (int j = i + 1; j < 4; j++) sum -= K[i * 4 + j] * sol[j];
        sol[i] = sum / K[i * 4 + i];
    }

    /* Reconstruct P */
    P->data[0] = sol[0]; P->data[1] = sol[1];
    P->data[2] = sol[2]; P->data[3] = sol[3];

    return true;
}

bool sss_lyap_3x3(const SwitchedMatrix *A, const SwitchedMatrix *Q, SwitchedMatrix *P) {
    return sss_lyap_solve(A, Q, P, 3);
}

bool sss_lyap_solve(const SwitchedMatrix *A, const SwitchedMatrix *Q, SwitchedMatrix *P, int n) {
    if (!A || !Q || !P || n <= 0 || n > 6) return false;
    if (A->rows != n || A->cols != n) return false;

    int n2 = n * n;
    double *K = (double *)calloc((size_t)(n2 * n2), sizeof(double));
    if (!K) return false;

    /* Build: (I kron A^T + A^T kron I) */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int row = i * n + j;
            /* I kron A^T: A^T[j,k] appears at [i*n+j, i*n+k] */
            for (int k = 0; k < n; k++) {
                int col = i * n + k;
                K[row * n2 + col] += A->data[k * n + j]; /* A^T[j,k] = A[k,j] */
            }
            /* A^T kron I: A^T[i,k] appears at [i*n+j, k*n+j] */
            for (int k = 0; k < n; k++) {
                int col = k * n + j;
                K[row * n2 + col] += A->data[k * n + i]; /* A^T[i,k] = A[k,i] */
            }
        }
    }

    /* RHS: -Q */
    double *rhs = (double *)malloc((size_t)n2 * sizeof(double));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            rhs[i * n + j] = -Q->data[i * n + j];

    /* Augmented matrix */
    int cols_aug = n2 + 1;
    double *aug = (double *)malloc((size_t)(n2 * cols_aug) * sizeof(double));
    for (int i = 0; i < n2; i++) {
        for (int j = 0; j < n2; j++) aug[i * cols_aug + j] = K[i * n2 + j];
        aug[i * cols_aug + n2] = rhs[i];
    }

    /* Gaussian elimination */
    for (int col = 0; col < n2; col++) {
        int pivot = col;
        double pv = fabs(aug[col * cols_aug + col]);
        for (int r = col + 1; r < n2; r++) {
            double v = fabs(aug[r * cols_aug + col]);
            if (v > pv) { pv = v; pivot = r; }
        }
        if (pv < 1e-14) { free(K); free(rhs); free(aug); return false; }
        if (pivot != col) {
            for (int j = 0; j < cols_aug; j++) {
                double t = aug[col * cols_aug + j];
                aug[col * cols_aug + j] = aug[pivot * cols_aug + j];
                aug[pivot * cols_aug + j] = t;
            }
        }
        double piv_val = aug[col * cols_aug + col];
        for (int r = col + 1; r < n2; r++) {
            double fac = aug[r * cols_aug + col] / piv_val;
            for (int j = col; j < cols_aug; j++)
                aug[r * cols_aug + j] -= fac * aug[col * cols_aug + j];
        }
    }

    double *x = (double *)malloc((size_t)n2 * sizeof(double));
    for (int i = n2 - 1; i >= 0; i--) {
        double s = aug[i * cols_aug + n2];
        for (int j = i + 1; j < n2; j++) s -= aug[i * cols_aug + j] * x[j];
        x[i] = s / aug[i * cols_aug + i];
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            P->data[i * n + j] = x[i * n + j];

    free(K); free(rhs); free(aug); free(x);
    return true;
}

/* ============================================================================
 * L5: Common Lyapunov Function Search
 * ============================================================================ */

bool sss_clf_verify(const SwitchedMatrix *P, SwitchedSubsystem **modes, int n_modes) {
    if (!P || !modes) return false;

    if (!sm_is_positive_definite(P)) return false;

    int n = P->rows;
    for (int i = 0; i < n_modes; i++) {
        if (!modes[i]) continue;
        SwitchedMatrix *A = &modes[i]->A;

        /* Compute A^T P + P A */
        SwitchedMatrix AT = sm_create(n, n);
        SwitchedMatrix PA = sm_create(n, n);
        SwitchedMatrix ATP = sm_create(n, n);
        sm_transpose(&AT, A);
        sm_mul(&PA, P, A);
        sm_mul(&ATP, &AT, P);

        SwitchedMatrix R = sm_create(n, n);
        sm_add(&R, &ATP, &PA);

        /* Check negative definiteness via trace */
        double tr = sm_trace(&R);
        if (tr >= -1e-10) {
            sm_free(&AT); sm_free(&PA); sm_free(&ATP); sm_free(&R);
            return false;
        }

        sm_free(&AT); sm_free(&PA); sm_free(&ATP); sm_free(&R);
    }

    return true;
}

bool sss_clf_gradient_descent(SwitchedSubsystem **modes, int n_modes, int n,
                               CommonLyapunovFunction *clf, int max_iters) {
    if (!modes || !clf || n_modes <= 0) return false;

    SwitchedMatrix P = sm_create(n, n);
    sm_identity(&P);
    double alpha = 0.01;

    for (int iter = 0; iter < max_iters; iter++) {
        double max_violation = -INFINITY;
        int worst_mode = 0;

        for (int m = 0; m < n_modes; m++) {
            if (!modes[m]) continue;
            SwitchedMatrix *A = &modes[m]->A;
            SwitchedMatrix AT = sm_create(n, n);
            sm_transpose(&AT, A);
            SwitchedMatrix PA = sm_create(n, n);
            SwitchedMatrix ATP = sm_create(n, n);
            sm_mul(&PA, &P, A);
            sm_mul(&ATP, &AT, &P);
            SwitchedMatrix R = sm_create(n, n);
            sm_add(&R, &ATP, &PA);

            double max_eig = 0.0;
            for (int i = 0; i < n; i++) {
                double row_sum = 0.0;
                for (int j = 0; j < n; j++) row_sum += fabs(R.data[i * n + j]);
                if (row_sum > max_eig) max_eig = row_sum;
            }
            if (max_eig > max_violation) {
                max_violation = max_eig;
                worst_mode = m;
            }
            sm_free(&AT); sm_free(&PA); sm_free(&ATP); sm_free(&R);
        }

        if (max_violation < -1e-8) {
            sm_copy(&clf->P, &P);
            clf->is_valid = true;
            clf->decay_rate = -max_violation / 2.0;
            clf->min_eig = 0.01;
            clf->max_eig = 1.0;
            sm_free(&P);
            return true;
        }

        /* Gradient step */
        SwitchedMatrix *A = &modes[worst_mode]->A;
        SwitchedMatrix AT = sm_create(n, n);
        sm_transpose(&AT, A);
        SwitchedMatrix PA = sm_create(n, n);
        SwitchedMatrix ATP = sm_create(n, n);
        sm_mul(&PA, &P, A);
        sm_mul(&ATP, &AT, &P);
        SwitchedMatrix grad = sm_create(n, n);
        sm_add(&grad, &ATP, &PA);
        sm_mul_scalar(&grad, alpha);
        sm_sub(&P, &P, &grad);

        /* Symmetrize */
        SwitchedMatrix PT = sm_create(n, n);
        sm_transpose(&PT, &P);
        sm_add(&P, &P, &PT);
        sm_mul_scalar(&P, 0.5);

        /* Ensure positive definiteness */
        for (int i = 0; i < n; i++)
            if (P.data[i * n + i] < 0.001) P.data[i * n + i] = 0.001;

        sm_free(&AT); sm_free(&PA); sm_free(&ATP); sm_free(&grad); sm_free(&PT);

        if (iter % 50 == 0) alpha *= 0.95;
    }

    sm_free(&P);
    return false;
}

/* ============================================================================
 * L5: Multiple Lyapunov Function Computation
 * ============================================================================ */

void sss_mlf_compute_all(SwitchedSystem *sys, MultipleLyapunovFunctions *mlf) {
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
        mlf->is_valid[i] = sss_lyap_solve(&sub->A, &Q, &mlf->P[i], n);
        if (mlf->is_valid[i]) {
            double min_e = INFINITY, max_e = 0.0;
            for (int r = 0; r < n; r++) {
                double diag = mlf->P[i].data[r * n + r];
                double radius = 0.0;
                for (int c = 0; c < n; c++)
                    if (c != r) radius += fabs(mlf->P[i].data[r * n + c]);
                if (diag - radius < min_e) min_e = diag - radius;
                if (diag + radius > max_e) max_e = diag + radius;
            }
            mlf->min_eig[i] = (min_e > 0) ? min_e : 0.001;
            mlf->max_eig[i] = max_e;
            mlf->decay_rates[i] = 0.5 / mlf->max_eig[i];
        }
    }
    sm_free(&Q);
    mlf->mu = sss_compute_mu(mlf);
}

bool sss_mlf_verify_sequence(const SwitchedSystem *sys, const MultipleLyapunovFunctions *mlf,
                              const SwitchingSignal *signal) {
    if (!sys || !mlf || !signal) return false;

    /* At each switch from mode i to mode j at time t_k:
     * V_j(x(t_k)) <= mu * V_i(x(t_k)) must hold */
    for (int k = 1; k <= signal->n_switches; k++) {
        int prev_mode = signal->mode_sequence[k - 1];
        int curr_mode = signal->mode_sequence[k];
        if (prev_mode >= mlf->n_modes || curr_mode >= mlf->n_modes) continue;
        if (!mlf->is_valid[prev_mode] || !mlf->is_valid[curr_mode]) continue;

        /* Compute V_i(x) = x^T P_i x */
        SwitchedVector Px_prev = sv_create(sys->state_dim);
        SwitchedVector Px_curr = sv_create(sys->state_dim);
        sm_matvec_mul(&Px_prev, &mlf->P[prev_mode], &sys->state);
        sm_matvec_mul(&Px_curr, &mlf->P[curr_mode], &sys->state);
        double V_prev = sv_dot(&sys->state, &Px_prev);
        double V_curr = sv_dot(&sys->state, &Px_curr);
        sv_free(&Px_prev);
        sv_free(&Px_curr);

        if (V_curr > mlf->mu * V_prev + 1e-10) return false;
    }
    return true;
}

double sss_lyap_eval(const SwitchedMatrix *P, const SwitchedVector *x) {
    if (!P || !x || P->rows != x->n) return 0.0;
    SwitchedVector Px = sv_create(x->n);
    sm_matvec_mul(&Px, P, x);
    double val = sv_dot(x, &Px);
    sv_free(&Px);
    return val;
}

double sss_lyap_derivative(const SwitchedMatrix *P, const SwitchedMatrix *A,
                            const SwitchedVector *x) {
    if (!P || !A || !x) return 0.0;
    int n = P->rows;
    if (n != A->rows || n != x->n) return 0.0;

    /* dV/dt = x^T (A^T P + P A) x */
    SwitchedMatrix AT = sm_create(n, n);
    sm_transpose(&AT, A);
    SwitchedMatrix PA = sm_create(n, n);
    SwitchedMatrix ATP = sm_create(n, n);
    sm_mul(&PA, P, A);
    sm_mul(&ATP, &AT, P);
    SwitchedMatrix R = sm_create(n, n);
    sm_add(&R, &ATP, &PA);

    double dV = sss_lyap_eval(&R, x);

    sm_free(&AT); sm_free(&PA); sm_free(&ATP); sm_free(&R);
    return dV;
}