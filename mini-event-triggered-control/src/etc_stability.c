#include "etc_stability.h"
#include "etc_trigger.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef ETC_ZENO_THRESHOLD
#define ETC_ZENO_THRESHOLD 1e6
#endif

/* ============================================================================
 * Comparison Functions (K, K∞, KL classes)
 * ============================================================================ */

double etc_kclass_function(double s, double a, double b, int form) {
    if (s < 0.0) s = 0.0;
    switch (form) {
        case 0: /* Linear: α(s) = a·s */
            return a * s;
        case 1: /* Quadratic: α(s) = a·s² */
            return a * s * s;
        case 2: /* Saturating: α(s) = a·s / (1 + b·s) */
            return a * s / (1.0 + b * s);
        default:
            return a * s;
    }
}

bool etc_verify_kclass(double (*alpha_fn)(double, double, double, int),
                        double a, double b, int form, double s_max) {
    if (!alpha_fn || a <= 0.0 || s_max <= 0.0) return false;

    /* Check α(0) = 0 */
    if (fabs(alpha_fn(0.0, a, b, form)) > 1e-12) return false;

    /* Check monotonicity: α(s1) < α(s2) for s1 < s2 */
    int n_samples = 100;
    double prev = alpha_fn(0.0, a, b, form);
    for (int i = 1; i <= n_samples; i++) {
        double s = (s_max * (double)i) / (double)n_samples;
        double curr = alpha_fn(s, a, b, form);
        if (curr < prev - 1e-12) return false;
        prev = curr;
    }
    return true;
}

double etc_kl_function(double s, double t, double k, double lam) {
    if (s < 0.0) s = 0.0;
    if (t < 0.0) t = 0.0;
    return k * s * exp(-lam * t);
}

/* ============================================================================
 * ISS Lyapunov Analysis
 *
 * Theorem (Sontag 2008, Tabuada 2007):
 *   Consider the ETC system:
 *     ẋ = Acl x − BK e
 *   with quadratic Lyapunov function V(x) = xᵀPx.
 *
 *   If P > 0 solves AclᵀP + PAcl = −Q with Q > 0, then:
 *     α₁(|x|) ≤ V(x) ≤ α₂(|x|)  with
 *       α₁(s) = λ_min(P)·s², α₂(s) = λ_max(P)·s²
 *
 *     V̇(x) = −xᵀQx + 2xᵀPBKe
 *          ≤ −λ_min(Q)|x|² + 2|PBK||x||e|
 *          ≤ −(1−θ)λ_min(Q)|x|² − θλ_min(Q)|x|² + 2|PBK||x||e|
 *
 *   Using Young's inequality: 2|PBK||x||e| ≤ (|PBK|²/θ)|x|² + θ|e|²
 *     V̇ ≤ −λ_min(Q)|x|² + (|PBK|²/θ)|x|² + θ|e|²
 *        ≤ −α₃(|x|) + γ(|e|)
 *   with α₃(s) = (λ_min(Q) − |PBK|²/θ)·s², γ(s) = θ·s².
 *
 *   ISS condition: λ_min(Q) > |PBK|²/θ for some θ > 0.
 *   By choosing θ = |PBK|, we get α₃(s) = (λ_min(Q) − |PBK|)·s².
 * ============================================================================ */

bool etc_verify_iss_lyapunov(const ETCSystem* sys,
                              double* alpha1, double* alpha2,
                              double* alpha3, double* gamma) {
    if (!sys || !alpha1 || !alpha2 || !alpha3 || !gamma) return false;

    /* For a linear system with quadratic V, the ISS-Lyapunov characterization
     * reduces to verifying the Lyapunov equation and computing gains. */
    const ETCMatrix* P = &sys->V.P;
    int n = sys->n_states;

    if (!P->data || n <= 0) return false;

    /* Compute P norm and estimate eigenvalues via Gershgorin */
    double diag_sum = 0.0;
    for (int i = 0; i < n; i++) diag_sum += P->data[i * n + i];
    double trace = diag_sum;
    double lambda_min_P = trace; /* Rough estimate — actual min may be smaller */
    double lambda_max_P = trace; /* Rough estimate — actual max may be larger */

    /* For a better estimate, use diagonal dominance */
    for (int i = 0; i < n; i++) {
        double radius = 0.0;
        for (int j = 0; j < n; j++)
            if (j != i) radius += fabs(P->data[i * n + j]);
        double li = P->data[i * n + i] - radius;
        double ui = P->data[i * n + i] + radius;
        if (li < lambda_min_P) lambda_min_P = li;
        if (ui > lambda_max_P) lambda_max_P = ui;
    }
    if (lambda_min_P < 1e-12) lambda_min_P = 1e-12;

    *alpha1 = lambda_min_P;
    *alpha2 = lambda_max_P;

    /* Compute ||PBK|| */
    ETCMatrix BK = etc_matrix_create(n, n);
    etc_matrix_mul(&sys->B, &sys->K, &BK);
    ETCMatrix PBK = etc_matrix_create(n, n);
    etc_matrix_mul(P, &BK, &PBK);
    double norm_PBK = etc_matrix_norm_2(&PBK);
    etc_matrix_free(&BK);
    etc_matrix_free(&PBK);

    /* Assume Q = I, so λ_min(Q) = 1.0 */
    double lambda_min_Q = 1.0;

    /* Check ISS condition: λ_min(Q) > ||PBK|| */
    if (lambda_min_Q <= norm_PBK) {
        *alpha3 = 0.0;
        *gamma = 0.0;
        return false;
    }

    *alpha3 = lambda_min_Q - norm_PBK;
    *gamma = norm_PBK; /* using θ = ||PBK|| */
    return true;
}

/* ============================================================================
 * Lyapunov Equation Solver
 *
 * Solves AᵀP + PA = −Q for P.
 * For small n (≤ 16): use vectorized Kronecker product approach.
 * For larger n: use gradient descent on the residual ||AᵀP + PA + Q||_F.
 * ============================================================================ */

bool etc_solve_lyapunov(const ETCMatrix* A, ETCMatrix* P,
                         const ETCMatrix* Q, int n) {
    if (!A || !P || !Q || n <= 0) return false;
    if (A->rows != n || A->cols != n) return false;

    /* For n = 2: analytic solution of A^T P + P A = -Q.
     * Let A = [a b; c d], Q = [q1 q2; q2 q3] (symmetric).
     * P = [p1 p2; p2 p3]. Then:
     *   2a·p1 + 2c·p2 = -q1
     *   b·p1 + (a+d)·p2 + c·p3 = -q2
     *   2b·p2 + 2d·p3 = -q3
     *
     * This is a 3×3 linear system:
     *   [2a   2c   0 ] [p1]   [-q1]
     *   [ b  a+d   c ] [p2] = [-q2]
     *   [ 0   2b  2d ] [p3]   [-q3]
     */
    if (n == 2) {
        double a = A->data[0], b = A->data[1];
        double c = A->data[2], d = A->data[3];
        double q1 = Q->data[0], q2 = Q->data[1], q3 = Q->data[3];

        /* Build 3x3 system M * p = r */
        double M[9] = {
            2.0*a, 2.0*c, 0.0,
            b, a + d, c,
            0.0, 2.0*b, 2.0*d
        };
        double r[3] = {-q1, -q2, -q3};

        /* Solve 3x3 via Cramer's rule */
        double detM = M[0]*(M[4]*M[8] - M[5]*M[7])
                    - M[1]*(M[3]*M[8] - M[5]*M[6])
                    + M[2]*(M[3]*M[7] - M[4]*M[6]);

        if (fabs(detM) < 1e-15) return false;

        double inv_det = 1.0 / detM;
        double p1 = inv_det * (r[0]*(M[4]*M[8] - M[5]*M[7])
                             - M[1]*(r[1]*M[8] - M[5]*r[2])
                             + M[2]*(r[1]*M[7] - M[4]*r[2]));
        double p2 = inv_det * (M[0]*(r[1]*M[8] - M[5]*r[2])
                             - r[0]*(M[3]*M[8] - M[5]*M[6])
                             + M[2]*(M[3]*r[2] - r[1]*M[6]));
        double p3 = inv_det * (M[0]*(M[4]*r[2] - r[1]*M[7])
                             - M[1]*(M[3]*r[2] - r[1]*M[6])
                             + r[0]*(M[3]*M[7] - M[4]*M[6]));

        P->data[0] = p1; P->data[1] = p2;
        P->data[2] = p2; P->data[3] = p3;
        return etc_matrix_is_positive_definite(P);
    }

    /* For n > 2: iterative gradient descent solver.
     * The continuous Lyapunov equation:
     *   Aᵀ P + P A = −Q
     *
     * P_{k+1} = P_k − η (AᵀP_k + P_kA + Q) */
    int n2 = n * n;
    double eta = 0.01;
    int max_iter = 5000;
    double tol = 1e-8;

    for (int i = 0; i < n2; i++) P->data[i] = Q->data[i];

    double* At = (double*)malloc((size_t)n2 * sizeof(double));
    double* AP = (double*)malloc((size_t)n2 * sizeof(double));
    double* PA = (double*)malloc((size_t)n2 * sizeof(double));
    double* residual = (double*)malloc((size_t)n2 * sizeof(double));
    if (!At || !AP || !PA || !residual) {
        free(At); free(AP); free(PA); free(residual);
        return false;
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            At[j * n + i] = A->data[i * n + j];

    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int k = 0; k < n; k++)
                    sum += A->data[i * n + k] * P->data[k * n + j];
                AP[i * n + j] = sum;
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int k = 0; k < n; k++)
                    sum += P->data[i * n + k] * A->data[k * n + j];
                PA[i * n + j] = sum;
            }
        }
        double frob_norm = 0.0;
        for (int i = 0; i < n2; i++) {
            int row = i / n, col = i % n;
            double atp = 0.0;
            for (int k = 0; k < n; k++)
                atp += At[row * n + k] * P->data[k * n + col];
            residual[i] = atp + PA[i] + Q->data[i];
            frob_norm += residual[i] * residual[i];
        }
        frob_norm = sqrt(frob_norm);
        if (frob_norm < tol) break;

        for (int i = 0; i < n2; i++)
            P->data[i] -= eta * residual[i];
    }

    free(At); free(AP); free(PA); free(residual);
    return etc_matrix_is_positive_definite(P);
}

bool etc_compute_clf(ETCSystem* sys) {
    if (!sys || sys->n_states <= 0) return false;
    int n = sys->n_states;

    /* Q = I */
    ETCMatrix Q = etc_matrix_create(n, n);
    for (int i = 0; i < n; i++) Q.data[i * n + i] = 1.0;

    bool ok = etc_solve_lyapunov(&sys->Acl, &sys->V.P, &Q, n);
    etc_matrix_free(&Q);

    if (ok) {
        sys->V.is_positive_definite = true;
        sys->V.P_norm = etc_matrix_norm_2(&sys->V.P);

        if (n == 2) {
            /* Exact 2x2 eigenvalues */
            double re[2], im[2];
            ETCMatrix Pcopy = etc_matrix_create(2, 2);
            for (int i = 0; i < 4; i++) Pcopy.data[i] = sys->V.P.data[i];
            etc_matrix_eigenvalues_2x2(&Pcopy, re, im);
            double lam1 = sqrt(re[0]*re[0] + im[0]*im[0]);
            double lam2 = sqrt(re[1]*re[1] + im[1]*im[1]);
            sys->V.lambda_min_P = (lam1 < lam2) ? lam1 : lam2;
            sys->V.lambda_max_P = (lam1 > lam2) ? lam1 : lam2;
            etc_matrix_free(&Pcopy);
        } else {
            /* Gershgorin estimate for larger matrices */
            double diag_sum = 0.0;
            for (int i = 0; i < n; i++) diag_sum += sys->V.P.data[i * n + i];
            sys->V.lambda_min_P = diag_sum > 0.0 ? diag_sum / (double)n : 0.0;
            sys->V.lambda_max_P = diag_sum;
            for (int i = 0; i < n; i++) {
                double radius = 0.0;
                for (int j = 0; j < n; j++)
                    if (j != i) radius += fabs(sys->V.P.data[i * n + j]);
                double li = sys->V.P.data[i * n + i] - radius;
                double ui = sys->V.P.data[i * n + i] + radius;
                if (li < sys->V.lambda_min_P) sys->V.lambda_min_P = li;
                if (ui > sys->V.lambda_max_P) sys->V.lambda_max_P = ui;
            }
            if (sys->V.lambda_min_P < 1e-12) sys->V.lambda_min_P = 1e-12;
        }
    }
    return ok;
}

/* ============================================================================
 * L2-Gain Analysis
 *
 * For the ETC system ẋ = Acl x − BK e + Bw w,
 * the L₂-gain from w to x satisfies:
 *
 *   ∫₀ᵀ |x|² dt ≤ γ² ∫₀ᵀ |w|² dt + β(|x(0)|)
 *
 * Under event-triggering with threshold σ, the error satisfies |e| ≤ σ|x|
 * between events. This induces an additional perturbation term in the
 * Lyapunov analysis.
 *
 * Using the ISS-Lyapunov function V(x) = xᵀPx:
 *   V̇ ≤ −|x|² + γ²|w|² + c|e|²
 *   with c = ||PBK||².
 *
 * Under the trigger condition |e| ≤ σ|x|:
 *   V̇ ≤ −(1 − c σ²)|x|² + γ²|w|²
 *
 * For stability we need 1 > c σ², i.e., σ < 1/√c.
 * The L₂-gain is then γ = ||PBw|| / (1 − c σ²).
 * ============================================================================ */

bool etc_compute_l2_gain(const ETCSystem* sys, const ETCMatrix* Bw,
                          double sigma, double* gamma) {
    if (!sys || !Bw || !gamma) return false;
    int n = sys->n_states;

    /* Compute ||PBK|| */
    ETCMatrix BK = etc_matrix_create(n, n);
    etc_matrix_mul(&sys->B, &sys->K, &BK);
    ETCMatrix PBK = etc_matrix_create(n, n);
    etc_matrix_mul(&sys->V.P, &BK, &PBK);
    double norm_PBK2 = etc_matrix_norm_2(&PBK);
    etc_matrix_free(&BK);
    etc_matrix_free(&PBK);
    norm_PBK2 = norm_PBK2 * norm_PBK2;

    /* Check stability condition */
    double c_sigma2 = norm_PBK2 * sigma * sigma;
    if (c_sigma2 >= 1.0) {
        *gamma = INFINITY;
        return false;
    }

    /* Compute ||PBw|| */
    ETCMatrix PBw = etc_matrix_create(n, Bw->cols);
    etc_matrix_mul(&sys->V.P, Bw, &PBw);
    double norm_PBw = etc_matrix_norm_2(&PBw);
    etc_matrix_free(&PBw);

    *gamma = norm_PBw / (1.0 - c_sigma2);
    return true;
}

/* ============================================================================
 * Zeno Analysis
 * ============================================================================ */

bool etc_check_zeno_free(const ETCSystem* sys) {
    if (!sys) return false;
    /* Theoretical check: for static triggering, Zeno-free iff σ < σ_max */
    ETCMatrix BK = etc_matrix_create(sys->n_states, sys->n_states);
    etc_matrix_mul(&sys->B, &sys->K, &BK);
    double sigma_max = etc_compute_sigma_max(&sys->V.P, &BK);
    etc_matrix_free(&BK);

    if (sys->sigma >= sigma_max) return false;
    if (sys->regime == ETC_REGIME_ZENO) return false;

    /* Empirical check from event history */
    if (sys->history.n_events > 10) {
        if (etc_history_check_zeno(&sys->history, ETC_ZENO_THRESHOLD))
            return false;
    }
    return true;
}

double etc_compute_zeno_time(const ETCSystem* sys) {
    if (!sys) return INFINITY;
    if (etc_check_zeno_free(sys)) return INFINITY;

    /* For systems that are not Zeno-free, estimate the Zeno time.
     * If σ ≥ σ_max, V̇ may be positive for some states, causing
     * accelerating event sequences.
     *
     * The Zeno time for linear ETC with σ > σ_max can be bounded by:
     * T_z ≤ 1 / (σ − σ_max) * (some constant)
     */
    ETCMatrix BK = etc_matrix_create(sys->n_states, sys->n_states);
    etc_matrix_mul(&sys->B, &sys->K, &BK);
    double sigma_max = etc_compute_sigma_max(&sys->V.P, &BK);
    etc_matrix_free(&BK);

    if (sys->sigma <= sigma_max) return INFINITY;

    double diff = sys->sigma - sigma_max;
    return 1.0 / (diff * 100.0); /* Rough estimate */
}

/* ============================================================================
 * Practical Stability
 * ============================================================================ */

bool etc_check_practical_stability(const ETCSystem* sys, double* radius) {
    if (!sys || !radius) return false;

    /* For ETC with |e| ≤ σ|x|, the system is practically stable if
     * trajectories converge to a ball of radius r around the origin.
     *
     * The ultimate bound is:
     *   r = σ ||BK|| λ_max(P) / (λ_min(Q) (1 − σ||BK||/λ_min(Q)))
     *
     * This follows from the Lyapunov analysis with the trigger condition.
     */
    int n = sys->n_states;

    /* Compute ||BK|| */
    ETCMatrix BK_mat = etc_matrix_create(n, n);
    etc_matrix_mul(&sys->B, &sys->K, &BK_mat);
    double norm_BK = etc_matrix_norm_2(&BK_mat);
    etc_matrix_free(&BK_mat);

    double sigma = sys->sigma;
    /* normalizer from ISS analysis */
    double alpha3 = 1.0; /* λ_min(Q) assumption */
    double lambda_max_P = sys->V.lambda_max_P;
    if (lambda_max_P < 1e-12) lambda_max_P = 1.0;

    /* Compute ultimate bound radius */
    double r = sigma * norm_BK * lambda_max_P;
    double denom = alpha3 * (1.0 - sigma * norm_BK / alpha3);
    if (denom > 1e-12)
        r = r / denom;
    else
        r = 1e6; /* Effectively unstable */

    *radius = r;

    /* Check if current state norm is within the ultimate bound */
    double x_norm = etc_system_state_norm(sys);
    return x_norm <= r * 1.1; /* 10% tolerance */
}

double etc_ultimate_bound(const ETCSystem* sys, double sigma) {
    if (!sys) return INFINITY;

    int n = sys->n_states;
    ETCMatrix BK_mat = etc_matrix_create(n, n);
    etc_matrix_mul(&sys->B, &sys->K, &BK_mat);
    double norm_BK = etc_matrix_norm_2(&BK_mat);
    etc_matrix_free(&BK_mat);

    double lambda_min_Q = 1.0; /* Q = I assumption */
    double lambda_max_P = sys->V.lambda_max_P;
    if (lambda_max_P < 1e-12) lambda_max_P = 1.0;

    double denom = (1.0 - sigma) * lambda_min_Q;
    if (denom <= 1e-15) return INFINITY;

    /* Ultimate bound: b = σ ||BK|| λ_max(P) / ((1−σ) λ_min(Q)) */
    return sigma * norm_BK * lambda_max_P / denom;
}
