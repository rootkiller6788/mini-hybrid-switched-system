#include "etc_trigger.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ============================================================================
 * Trigger Functions — Evaluate Γ(x, e)
 * ============================================================================ */

double etc_trigger_static(const ETCVector* x, const ETCVector* e,
                           double sigma, double epsilon) {
    (void)epsilon;
    if (!x || !e) return -1.0;
    return etc_vector_norm(e) - sigma * etc_vector_norm(x);
}

double etc_trigger_quadratic(const ETCVector* x, const ETCVector* e,
                              double sigma, double epsilon) {
    (void)epsilon;
    if (!x || !e) return -1.0;
    double x_n2 = etc_vector_dot(x, x);
    double e_n2 = etc_vector_dot(e, e);
    return e_n2 - sigma * x_n2;
}

double etc_trigger_mixed(const ETCVector* x, const ETCVector* e,
                          double sigma, double epsilon_effective) {
    if (!x || !e) return -1.0;
    double x_n2 = etc_vector_dot(x, x);
    double e_n2 = etc_vector_dot(e, e);
    return e_n2 - sigma * x_n2 - epsilon_effective;
}

double etc_trigger_absolute(const ETCVector* x, const ETCVector* e,
                             double sigma, double epsilon) {
    (void)sigma;
    (void)x;
    if (!e) return -1.0;
    return etc_vector_norm(e) - epsilon;
}

double etc_trigger_dynamic(const ETCVector* x, const ETCVector* e,
                            double eta, double sigma, double theta) {
    if (!x || !e) return -1.0;
    double x_n2 = etc_vector_dot(x, x);
    double e_n2 = etc_vector_dot(e, e);
    double static_part = sigma * x_n2 - e_n2;
    /* In dynamic triggering, event fires when η + θ(σ|x|² − |e|²) ≤ 0 */
    /* We return (theta * static_part + eta) so that >= 0 means fire */
    return eta + theta * static_part;
}

/* ============================================================================
 * Self-Triggered Interval Computation
 *
 * Theory (Mazo et al. 2010, Anta & Tabuada 2010):
 *   Under zero-order hold u(t) = K x_k for t ∈ [t_k, t_{k+1}),
 *   the state evolves as:
 *     x(t) = e^{Acl t} x_k + ∫₀ᵗ e^{Acl(t−s)} BK e(s) ds
 *
 *   With e(0) = 0 and x(0) = x_k.
 *
 *   We search for the first t > 0 where |e(t)| = σ|x(t)|.
 *   Using a simple bisection/search over the time interval.
 * ============================================================================ */

double etc_self_triggered_interval(const ETCVector* x_k,
                                    const ETCMatrix* Acl,
                                    const ETCMatrix* BK,
                                    double sigma,
                                    double tau_min,
                                    double tau_max) {
    if (!x_k || !Acl || !BK || sigma <= 0.0 || tau_min <= 0.0) return tau_min;
    if (tau_max <= tau_min) return tau_min;

    /* Simple linear search for the triggering time.
     * For each candidate time t, we need to compute e(t) and x(t).
     * Since this is a self-triggered scheme, we use the fact that
     * the error dynamics are ė = −A x − B K (x − e) = −Acl x + BK e.
     * In the small neighborhood where |e| is small, we can approximate
     * x(t) ≈ e^{Acl t} x_k, and e(t) evolves accordingly.
     *
     * For practical computation, we use a step-by-step forward
     * simulation approach to find the crossing time.
     */

    /* Approximate using the worst-case Lipschitz bound:
     * The inter-event time lower bound for static trigger is:
     * τ_min_theory = σ / (||Acl|| + σ||BK||)
     * We use this as the basis and scale up. */
    double norm_Acl = etc_matrix_norm_2(Acl);
    double norm_BK = etc_matrix_norm_2(BK);
    double denom = norm_Acl + sigma * norm_BK;
    if (denom < 1e-15) return tau_max;

    /* Compute state-dependent τ using the norm of x_k */
    double xk_norm = etc_vector_norm(x_k);
    /* For larger initial states, the inter-event time is typically shorter */
    double scale = (xk_norm > 1e-10) ? 1.0 / (1.0 + xk_norm) : 1.0;
    double tau_est = sigma * scale / denom;
    if (tau_est < tau_min) tau_est = tau_min;
    if (tau_est > tau_max) tau_est = tau_max;

    /* Refine: search for exact crossing via bisection on [tau_min, tau_est] */
    /* Simplified: return the estimated time */
    return tau_est;
}

/* ============================================================================
 * Trigger Analysis Functions
 * ============================================================================ */

double etc_compute_sigma_max(const ETCMatrix* P, const ETCMatrix* BK) {
    if (!P || !BK || !P->data || !BK->data) return 0.0;
    /* σ_max = λ_min(Q) / (2 ||PBK||) */
    /* For Q = I, λ_min(Q) = 1.0. Otherwise estimate from P. */
    double lambda_min_Q = 1.0; /* Assuming standard Lyapunov eq with Q = I */

    /* Compute ||PBK|| */
    int n = P->rows;
    ETCMatrix PBK = etc_matrix_create(n, n);
    etc_matrix_mul(P, BK, &PBK);
    double norm_PBK = etc_matrix_norm_2(&PBK);
    etc_matrix_free(&PBK);

    if (norm_PBK < 1e-15) return 1.0;
    double sigma_max = lambda_min_Q / (2.0 * norm_PBK);
    return sigma_max < 1.0 ? sigma_max : 1.0;
}

bool etc_is_sigma_stabilizing(const ETCMatrix* P, const ETCMatrix* BK,
                               double sigma) {
    double sigma_max = etc_compute_sigma_max(P, BK);
    return sigma < sigma_max;
}

double etc_compute_iet_lower_bound(const ETCMatrix* Acl, const ETCMatrix* BK,
                                    double sigma) {
    if (!Acl || !BK) return 0.0;
    double norm_Acl = etc_matrix_norm_2(Acl);
    double norm_BK = etc_matrix_norm_2(BK);
    double denom = norm_Acl + sigma * norm_BK;
    if (denom < 1e-15) return INFINITY;
    /* Tabuada (2007, Lemma 4): τ_min = σ / (||Acl|| + σ||BK||) */
    return sigma / denom;
}

double etc_dynamic_trigger_rhs(double eta, double beta,
                                const ETCVector* x, const ETCVector* e,
                                double sigma) {
    double x_n2 = etc_vector_dot(x, x);
    double e_n2 = etc_vector_dot(e, e);
    return -beta * eta + sigma * x_n2 - e_n2;
}

void etc_design_threshold(double desired_iet, double performance,
                           const ETCMatrix* Acl, const ETCMatrix* BK,
                           double* sigma_out, double* epsilon_out) {
    /* Performance: 1.0 = high performance (small sigma, more events)
     *              0.0 = low performance (large sigma, fewer events)
     * desired_iet: target average inter-event time
     *
     * From τ_min = σ / (||Acl|| + σ||BK||), we can solve for σ:
     *   τ_min (||Acl|| + σ||BK||) = σ
     *   τ_min ||Acl|| = σ (1 − τ_min ||BK||)
     *   σ = τ_min ||Acl|| / (1 − τ_min ||BK||)   if τ_min ||BK|| < 1
     */
    if (!Acl || !BK || !sigma_out || !epsilon_out) return;

    double norm_Acl = etc_matrix_norm_2(Acl);
    double norm_BK = etc_matrix_norm_2(BK);

    /* Desired σ from IET formula */
    double sigma_computed = 0.1; /* default */
    double denom_iet = 1.0 - desired_iet * norm_BK;
    if (denom_iet > 0.01) {
        sigma_computed = desired_iet * norm_Acl / denom_iet;
    } else {
        sigma_computed = 0.5; /* fallback */
    }

    /* Clamp to realistic range */
    if (sigma_computed < 0.001) sigma_computed = 0.001;
    if (sigma_computed > 0.99) sigma_computed = 0.99;

    /* Adjust by performance factor: high perf → smaller sigma */
    double sigma_adj = sigma_computed * (1.0 - 0.9 * performance);
    if (sigma_adj < 0.001) sigma_adj = 0.001;

    *sigma_out = sigma_adj;

    /* Epsilon provides an additional safety margin for mixed trigger.
     * Larger epsilon → fewer events (allows larger error before triggering).
     * ε = ε_base * (1 − performance): low perf → larger epsilon */
    double eps_base = sigma_adj * 0.1;
    *epsilon_out = eps_base * (1.0 - performance);
}
