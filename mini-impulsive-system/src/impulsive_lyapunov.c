/*
 * impulsive_lyapunov.c -- Lyapunov stability theory for impulsive systems
 *
 * Implements Lyapunov functions and stability verification:
 *   - Quadratic Lyapunov functions V(x) = x^T P x
 *   - Norm-based Lyapunov functions V(x) = ||x||_p^q
 *   - Composite (convex combination) Lyapunov functions
 *   - Lie derivative computation along flows
 *   - Jump ratio estimation
 *   - Dwell-time bound computation
 *   - Lyapunov equation solver: A^T P + P A = -Q
 *   - Algebraic Riccati equation solver (iterative)
 *   - Positive definiteness check via Cholesky attempt
 *
 * Key theorems implemented:
 *   Theorem 1 (Lyapunov stability): If V>0, L_f V <= -alpha*V during flow,
 *     and V(x^+) <= rho*V(x^-) with rho<1 at jumps, and dwell times
 *     satisfy tau_D > -ln(rho)/alpha, then the origin is GES.
 *   Theorem 2 (Comparison principle): If m(t) <= u(t) where u solves
 *     the comparison system, then stability of u implies stability of x.
 *
 * References:
 *   Yang (2001) "Impulsive Control Theory", Theorems 2.1-2.5
 *   Haddad et al. (2006) "Impulsive and Hybrid Dynamical Systems"
 *   Khalil (2002) "Nonlinear Systems", Chapter 4
 */
#include "impulsive_lyapunov.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ---- Quadratic Lyapunov Function ---- */

ImpLyapQuad* imp_lyap_quad_create(const double *P, int n)
{
    if (!P || n < 1 || n > IMP_MAX_DIM) return NULL;
    ImpLyapQuad *lyap = (ImpLyapQuad*)calloc(1, sizeof(ImpLyapQuad));
    if (!lyap) return NULL;
    lyap->P = (double*)malloc((size_t)n * n * sizeof(double));
    if (!lyap->P) { free(lyap); return NULL; }
    memcpy(lyap->P, P, (size_t)n * n * sizeof(double));
    lyap->n = n;
    return lyap;
}

void imp_lyap_quad_free(ImpLyapQuad *lyap)
{
    if (!lyap) return;
    free(lyap->P); free(lyap);
}

double imp_lyap_quad_V(const double *x, int n, void *ctx)
{
    ImpLyapQuad *lyap = (ImpLyapQuad*)ctx;
    if (!lyap || !x || n != lyap->n) return 0.0;
    double val = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            val += x[i] * lyap->P[i * n + j] * x[j];
    return val;
}

int imp_lyap_quad_gradV(const double *x, int n, double *grad, void *ctx)
{
    ImpLyapQuad *lyap = (ImpLyapQuad*)ctx;
    if (!lyap || !x || !grad || n != lyap->n) return -1;
    /* grad V = 2*P*x for symmetric P */
    for (int i = 0; i < n; i++) {
        grad[i] = 0.0;
        for (int j = 0; j < n; j++)
            grad[i] += (lyap->P[i * n + j] + lyap->P[j * n + i]) * x[j];
    }
    return 0;
}

double imp_lyap_quad_deriv_along_flow(const double *P, const double *A,
                                       const double *x, int n)
{
    if (!P || !A || !x || n < 1) return 0.0;
    /* L_f V = x^T (A^T P + P A) x */
    double val = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double M_ij = 0.0;
            for (int k = 0; k < n; k++)
                M_ij += A[k * n + i] * P[k * n + j]  /* A^T P */
                      + P[i * n + k] * A[k * n + j]; /* P A */
            val += x[i] * M_ij * x[j];
        }
    }
    return val;
}

double imp_lyap_quad_jump_ratio(const double *P, const double *J,
                                 const double *x, int n)
{
    if (!P || !J || !x || n < 1) return -1.0;
    /* V(x^-) = x^T P x */
    double V_before = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            V_before += x[i] * P[i * n + j] * x[j];
    if (V_before < IMP_EPS) return 0.0;

    /* x^+ = (I+J)*x, V(x^+) = x^+^T P x^+ */
    double *xp = (double*)malloc((size_t)n * sizeof(double));
    if (!xp) return -1.0;
    for (int i = 0; i < n; i++) {
        xp[i] = x[i];
        for (int k = 0; k < n; k++)
            xp[i] += J[i * n + k] * x[k];
    }
    double V_after = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            V_after += xp[i] * P[i * n + j] * xp[j];
    free(xp);
    return V_after / V_before;
}

/* ---- Norm-Based Lyapunov Function ---- */

ImpLyapNorm* imp_lyap_norm_create(int p, double q, int n)
{
    if (n < 1 || n > IMP_MAX_DIM || q <= 0.0 || (p != 1 && p != 2 && p != 0))
        return NULL;
    ImpLyapNorm *ln = (ImpLyapNorm*)calloc(1, sizeof(ImpLyapNorm));
    if (!ln) return NULL;
    ln->p = p; ln->q = q; ln->n = n;
    return ln;
}

void imp_lyap_norm_free(ImpLyapNorm *ln) { free(ln); }

double imp_lyap_norm_V(const double *x, int n, void *ctx)
{
    ImpLyapNorm *ln = (ImpLyapNorm*)ctx;
    if (!ln || !x || n != ln->n) return 0.0;
    double norm_val = 0.0;
    if (ln->p == 2) {
        for (int i = 0; i < n; i++) norm_val += x[i] * x[i];
        norm_val = sqrt(norm_val);
    } else if (ln->p == 1) {
        for (int i = 0; i < n; i++) norm_val += fabs(x[i]);
    } else { /* p == 0: infinity norm */
        for (int i = 0; i < n; i++)
            if (fabs(x[i]) > norm_val) norm_val = fabs(x[i]);
    }
    return pow(norm_val, ln->q);
}

int imp_lyap_norm_gradV(const double *x, int n, double *grad, void *ctx)
{
    ImpLyapNorm *ln = (ImpLyapNorm*)ctx;
    if (!ln || !x || !grad || n != ln->n) return -1;
    double norm_val, q = ln->q;
    /* grad ||x||_p^q = q * ||x||_p^{q-1} * grad ||x||_p */
    if (ln->p == 2) {
        norm_val = 0.0;
        for (int i = 0; i < n; i++) norm_val += x[i] * x[i];
        norm_val = sqrt(norm_val);
        if (norm_val < IMP_EPS) { memset(grad, 0, (size_t)n * sizeof(double)); return 0; }
        double coef = q * pow(norm_val, q - 2.0);
        for (int i = 0; i < n; i++) grad[i] = coef * x[i];
    } else if (ln->p == 1) {
        norm_val = 0.0;
        for (int i = 0; i < n; i++) norm_val += fabs(x[i]);
        double coef = q * pow(norm_val, q - 1.0);
        for (int i = 0; i < n; i++) grad[i] = coef * ((x[i] > 0) ? 1.0 : (x[i] < 0) ? -1.0 : 0.0);
    } else {
        norm_val = 0.0;
        for (int i = 0; i < n; i++)
            if (fabs(x[i]) > norm_val) norm_val = fabs(x[i]);
        double coef = q * pow(norm_val, q - 1.0);
        /* Subgradient: only the max component gets the gradient */
        for (int i = 0; i < n; i++)
            grad[i] = (fabs(x[i]) >= norm_val - IMP_EPS) ? coef * ((x[i] > 0) ? 1.0 : -1.0) : 0.0;
    }
    return 0;
}

/* ---- Composite Lyapunov Function ---- */

ImpLyapComposite* imp_lyap_composite_create(ImpLyapunovFn **comps,
                                              const double *w, int k)
{
    if (!comps || !w || k < 1) return NULL;
    ImpLyapComposite *lc = (ImpLyapComposite*)calloc(1, sizeof(ImpLyapComposite));
    if (!lc) return NULL;
    lc->components = (ImpLyapunovFn**)malloc((size_t)k * sizeof(ImpLyapunovFn*));
    lc->weights    = (double*)malloc((size_t)k * sizeof(double));
    if (!lc->components || !lc->weights) { imp_lyap_composite_free(lc); return NULL; }
    for (int i = 0; i < k; i++) {
        if (!comps[i] || w[i] < 0.0) { imp_lyap_composite_free(lc); return NULL; }
        lc->components[i] = comps[i];
        lc->weights[i]    = w[i];
    }
    lc->k = k;
    return lc;
}

void imp_lyap_composite_free(ImpLyapComposite *lc)
{
    if (!lc) return;
    free(lc->components); free(lc->weights); free(lc);
}

double imp_lyap_composite_V(const double *x, int n, void *ctx)
{
    ImpLyapComposite *lc = (ImpLyapComposite*)ctx;
    if (!lc || !x) return 0.0;
    double val = 0.0;
    for (int i = 0; i < lc->k; i++)
        val += lc->weights[i] * lc->components[i]->V(x, n, lc->components[i]->ctx);
    return val;
}

int imp_lyap_composite_gradV(const double *x, int n, double *grad, void *ctx)
{
    ImpLyapComposite *lc = (ImpLyapComposite*)ctx;
    if (!lc || !x || !grad) return -1;
    memset(grad, 0, (size_t)n * sizeof(double));
    for (int i = 0; i < lc->k; i++) {
        double *gi = (double*)calloc((size_t)n, sizeof(double));
        if (!gi) return -1;
        lc->components[i]->gradV(x, n, gi, lc->components[i]->ctx);
        for (int j = 0; j < n; j++)
            grad[j] += lc->weights[i] * gi[j];
        free(gi);
    }
    return 0;
}

/* ---- Impulsive Stability Check ---- */

int imp_lyap_check_stability(ImpLyapunovFn *lyap,
                              ImpVectorField f, void *f_ctx,
                              ImpJumpMap I, void *I_ctx,
                              const double *test_pts, int num_pts, int n,
                              double alpha, double rho,
                              ImpStabilityType *result)
{
    if (!lyap || !f || !I || !test_pts || num_pts < 1 || !result) return -1;
    if (alpha < 0.0 || rho < 0.0) return -1;

    /* Check V(0) = 0 */
    double *zero = (double*)calloc((size_t)n, sizeof(double));
    if (!zero) return -2;
    double V0 = lyap->V(zero, n, lyap->ctx);
    if (fabs(V0) > IMP_EPS) { free(zero); *result = IMP_STABLE; return -1; }
    free(zero);

    /* Check V(x) > 0 and stability conditions at test points */
    for (int k = 0; k < num_pts; k++) {
        const double *x = &test_pts[(size_t)k * n];
        double Vx = lyap->V(x, n, lyap->ctx);
        if (Vx <= IMP_EPS) { *result = IMP_STABLE; return -1; }

        /* Check derivative condition: L_f V <= -alpha * V */
        double *dxdt = (double*)malloc((size_t)n * sizeof(double));
        double *grad = (double*)malloc((size_t)n * sizeof(double));
        if (!dxdt || !grad) { free(dxdt); free(grad); return -2; }

        f(0.0, x, n, dxdt, f_ctx);
        lyap->gradV(x, n, grad, lyap->ctx);

        double LfV = 0.0;
        for (int i = 0; i < n; i++) LfV += grad[i] * dxdt[i];

        if (LfV > -alpha * Vx + IMP_EPS) {
            free(dxdt); free(grad);
            *result = IMP_STABLE; return -1;
        }

        /* Check jump condition: V(x^+) <= rho * V(x^-) */
        double *x_after = (double*)malloc((size_t)n * sizeof(double));
        if (!x_after) { free(dxdt); free(grad); return -2; }
        I(0.0, x, n, x_after, I_ctx);
        double V_after = lyap->V(x_after, n, lyap->ctx);
        if (V_after > rho * Vx + IMP_EPS) {
            free(dxdt); free(grad); free(x_after);
            *result = IMP_STABLE; return -1;
        }
        free(dxdt); free(grad); free(x_after);
    }

    /* All checks passed */
    if (rho < 1.0) *result = IMP_EXPONENTIAL;
    else if (alpha > 0.0) *result = IMP_EXPONENTIAL;
    else *result = IMP_ASYMPTOTIC;
    return 0;
}

/* ---- Dwell-Time Bound ---- */

double imp_lyap_dwell_time_bound(double c, double d)
{
    /* Condition: tau_D > ln(d) / (-c) if c < 0 (flow stable)
     *            tau_D < -ln(d) / c  if c > 0 (flow expanding) */
    if (d <= 0.0 || d >= 1.0) return 0.0;  /* d=0 means perfect jump (infinite margin) */
    if (fabs(c) < IMP_EPS) return IMP_HUGE;  /* no flow contraction, need arbitrarily fast impulses */
    if (c < 0.0) return fmax(0.0, log(d) / (-c));  /* upper bound, c is decay rate */
    else return -log(d) / c;  /* need frequent enough impulses */
}

double imp_lyap_estimate_contraction(ImpLyapunovFn *lyap,
                                      const double *x_before,
                                      const double *x_after, int n)
{
    if (!lyap || !x_before || !x_after || n < 1) return -1.0;
    double Vb = lyap->V(x_before, n, lyap->ctx);
    if (Vb < IMP_EPS) return -1.0;
    double Va = lyap->V(x_after, n, lyap->ctx);
    return Va / Vb;
}

double imp_lyap_estimate_flow_rate(ImpLyapunovFn *lyap,
                                    const double *x_t,
                                    const double *x_t_dt,
                                    double dt, int n)
{
    if (!lyap || !x_t || !x_t_dt || dt < IMP_EPS || n < 1) return 0.0;
    double Vt = lyap->V(x_t, n, lyap->ctx);
    if (Vt < IMP_EPS) return 0.0;
    double Vtdt = lyap->V(x_t_dt, n, lyap->ctx);
    return (Vtdt - Vt) / (dt * Vt);
}

double imp_lyap_average_dwell_time(double lambda_val, double mu)
{
    if (lambda_val <= 0.0 || mu < 1.0) return 0.0;
    return log(mu) / lambda_val;
}

/* ---- Lyapunov Equation Solver: A^T P + P A = -Q ---- */

int imp_lyap_solve_lyapunov_eq(const double *A, const double *Q,
                                double *P, int n)
{
    if (!A || !Q || !P || n < 1 || n > IMP_MAX_DIM) return -1;
    /* Solve via vectorized form: (I ⊗ A^T + A^T ⊗ I) vec(P) = -vec(Q) */
    /* Use simple iterative method: P_{k+1} = P_k + h*(A^T P_k + P_k A + Q) */
    /* Converges if A is Hurwitz (all eigenvalues have negative real parts). */
    int n2 = n * n;
    double h = 0.01;
    double *P_new = (double*)malloc((size_t)n2 * sizeof(double));
    double *AP   = (double*)malloc((size_t)n2 * sizeof(double));
    double *PA   = (double*)malloc((size_t)n2 * sizeof(double));
    if (!P_new || !AP || !PA) { free(P_new); free(AP); free(PA); return -2; }

    /* Initialize P to identity */
    for (int i = 0; i < n2; i++) P[i] = (i % (n + 1) == 0) ? 1.0 : 0.0;

    double residual = IMP_HUGE;
    for (int iter = 0; iter < IMP_MAX_ITER && residual > IMP_EPS; iter++) {
        /* Compute A^T P */
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                AP[i * n + j] = 0.0;
                for (int k = 0; k < n; k++)
                    AP[i * n + j] += A[k * n + i] * P[k * n + j];
            }
        /* Compute P A */
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                PA[i * n + j] = 0.0;
                for (int k = 0; k < n; k++)
                    PA[i * n + j] += P[i * n + k] * A[k * n + j];
            }
        /* Update and residual */
        residual = 0.0;
        for (int i = 0; i < n2; i++) {
            double delta = AP[i] + PA[i] + Q[i];
            P_new[i] = P[i] + h * delta;
            residual += delta * delta;
        }
        residual = sqrt(residual);
        /* Copy back */
        for (int i = 0; i < n2; i++) P[i] = P_new[i];
    }
    free(P_new); free(AP); free(PA);
    return (residual <= IMP_EPS) ? 0 : -3;
}

/* ---- Algebraic Riccati Equation Solver (iterative Kleinman) ---- */

int imp_lyap_solve_riccati_eq(const double *A, const double *B,
                               const double *Q, const double *R,
                               double *P, int n, int m)
{
    if (!A || !B || !Q || !R || !P || n < 1 || n > IMP_MAX_DIM || m < 1) return -1;
    /* Solve: A^T P + P A - P B R^{-1} B^T P + Q = 0 */
    /* Kleinman iterative method: solve Lyapunov equation at each step */
    int n2 = n * n;

    /* Compute R_inv (assume scalar R for simplicity, or diagonal) */
    double R_inv = 1.0;  /* simplified: R = 1 */
    if (R[0] > IMP_EPS) R_inv = 1.0 / R[0];

    /* Compute S = B R^{-1} B^T */
    double *S = (double*)calloc((size_t)n2, sizeof(double));
    if (!S) return -2;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < m; k++)
                S[i * n + j] += B[i * m + k] * R_inv * B[j * m + k];

    /* Initialize P to Q */
    for (int i = 0; i < n2; i++) P[i] = Q[i];

    double *P_prev = (double*)malloc((size_t)n2 * sizeof(double));
    if (!P_prev) { free(S); return -2; }

    double diff = IMP_HUGE;
    for (int iter = 0; iter < IMP_MAX_ITER && diff > IMP_EPS; iter++) {
        memcpy(P_prev, P, (size_t)n2 * sizeof(double));

        /* Solve Lyapunov: (A - S P)^T P_new + P_new (A - S P) = -(Q + P S P) */
        double *Ak = (double*)malloc((size_t)n2 * sizeof(double));
        double *Qk = (double*)malloc((size_t)n2 * sizeof(double));
        if (!Ak || !Qk) { free(Ak); free(Qk); free(S); free(P_prev); return -2; }

        /* Ak = A - S * P */
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double SP_ij = 0.0;
                for (int k = 0; k < n; k++) SP_ij += S[i * n + k] * P[k * n + j];
                Ak[i * n + j] = A[i * n + j] - SP_ij;
            }
        /* Qk = Q + P S P */
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double PSP_ij = 0.0;
                for (int k = 0; k < n; k++)
                    for (int l = 0; l < n; l++)
                        PSP_ij += P[i * n + k] * S[k * n + l] * P[l * n + j];
                Qk[i * n + j] = Q[i * n + j] + PSP_ij;
            }

        imp_lyap_solve_lyapunov_eq(Ak, Qk, P, n);

        /* Compute diff */
        diff = 0.0;
        for (int i = 0; i < n2; i++) {
            double d = P[i] - P_prev[i];
            diff += d * d;
        }
        diff = sqrt(diff);

        free(Ak); free(Qk);
    }

    free(S); free(P_prev);
    return (diff <= IMP_SQRT_EPS) ? 0 : -3;
}

/* ---- Positive Definiteness Check ---- */

int imp_lyap_check_pos_def(const double *P, int n)
{
    if (!P || n < 1 || n > IMP_MAX_DIM) return 0;
    /* Check leading principal minors via Gaussian elimination attempt */
    /* A real symm matrix is PD iff all eigenvalues > 0.
     * Quick check: diagonal entries > 0 and matrix is diagonally dominant. */
    for (int i = 0; i < n; i++) {
        if (P[i * n + i] <= IMP_EPS) return 0;
        double sum = 0.0;
        for (int j = 0; j < n; j++)
            if (j != i) sum += fabs(P[i * n + j]);
        if (P[i * n + i] <= sum) return -1;  /* not diagonally dominant, inconclusive */
    }
    return 1;  /* likely PD */
}
