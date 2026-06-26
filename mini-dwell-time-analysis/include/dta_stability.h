#ifndef DTA_STABILITY_H
#define DTA_STABILITY_H
#include "dta_core.h"
#include "dta_switch_signal.h"

/* ==============================================================
 * dta_stability.h - Stability Analysis of Switched Systems
 *
 * Core theorem: For a switched linear system ẋ = A_σ x, if there
 * exists a set of positive definite matrices {P_1, ..., P_m} and
 * scalars λ_i > 0, μ ≥ 1 such that:
 *
 *   A_i^T P_i + P_i A_i + 2λ_i P_i < 0   ∀i
 *   P_i ≤ μ P_j                            ∀i,j
 *
 * Then for any dwell time τ_d > ln(μ)/(2 min λ_i), the switched
 * system is GUES.
 *
 * Average dwell time variant (Hespanha & Morse 1999):
 *   τ_a > τ_a* = ln(μ) / λ,  λ = min λ_i
 *
 * References:
 *   Liberzon (2003) "Switching in Systems and Control", Ch. 3
 *   Hespanha & Morse (1999) IEEE CDC
 *   Geromel & Colaneri (2006) "Stability and stabilization of
 *     continuous-time switched linear systems", SIAM JCO
 * ============================================================== */

/* --- Dwell-time stability certificate --- */
typedef struct {
    DTA_StabilityVerdict verdict;
    double* lambda_i;          /* Decay rates per mode, length n_modes */
    double mu;                 /* Coupling constant: P_i ≤ μ P_j */
    double tau_d_star;         /* Theoretical minimum dwell time */
    double tau_d_actual;       /* Actual dwell time used */
    double decay_rate;         /* Guaranteed convergence rate */
    bool common_lyapunov;      /* True if common V exists */
    double* common_P;          /* Common Lyapunov matrix if exists, n×n */
    int state_dim;
} DTA_DwellStabilityResult;

/* --- Stability under arbitrary switching --- */
typedef struct {
    bool is_stable;
    bool common_lyapunov_exists;
    double* common_P;
    double min_eig;
    double max_eig;
    int n;
} DTA_ArbitrarySwitchingResult;

/* --- API: Stability Analysis --- */

/** Analyze GAS under given dwell time τ_d.
 *  Uses the multiple Lyapunov function approach.
 *  Returns a DTA_DwellStabilityResult with verdict and certificate. */
DTA_DwellStabilityResult dta_analyze_dwell_stability(
    const DTA_SwitchedSystem* sys, double tau_d);

/** Analyze stability under average dwell time τ_a with chatter N0.
 *  Implements Hespanha & Morse (1999) Theorem 2. */
DTA_DwellStabilityResult dta_analyze_adt_stability(
    const DTA_SwitchedSystem* sys, double tau_a, double N0);

/** Check if a switched linear system is stable under arbitrary switching.
 *  This requires a common Lyapunov function V(x)=x^T P x with
 *  A_i^T P + P A_i < 0 for all i.
 *  Solves via convex combination / iterative LMI approach. */
DTA_ArbitrarySwitchingResult dta_check_arbitrary_switching(
    const DTA_SwitchedSystem* sys);

/** Compute the minimum dwell time τ_d* that guarantees GAS.
 *  Uses bisection search with LMI feasibility checks.
 *  Returns INFINITY if no finite τ_d* can guarantee stability. */
double dta_compute_min_dwell_time(const DTA_SwitchedSystem* sys,
                                   double tau_max, double tol);

/** Compute the minimum average dwell time τ_a* that guarantees GAS. */
double dta_compute_min_avg_dwell_time(const DTA_SwitchedSystem* sys,
                                       double N0, double tau_max, double tol);

/** Verify the multiple Lyapunov function condition:
 *  P_i > 0, A_i^T P_i + P_i A_i + 2λ_i P_i < 0, P_i ≤ μ P_j */
bool dta_verify_mlf_condition(const DTA_SwitchedSystem* sys,
    const double** P_matrices, const double* lambda, double mu,
    int n_modes, int n);

/** Stability margin: largest γ such that A_i + γ I are all Hurwitz */
double dta_stability_margin(const DTA_SwitchedSystem* sys);

/** Exponential decay rate bound: ||x(t)|| ≤ c e^{-α t} ||x(0)|| */
double dta_decay_rate_bound(const DTA_SwitchedSystem* sys,
                             double tau_d);

/** Check if the switched system admits a common quadratic Lyapunov function.
 *  Uses the iterative algorithm from Liberzon (2003) §2.3. */
bool dta_find_common_quadratic(const DTA_SwitchedSystem* sys,
                                double* P_out);

/** Stability of switched nonlinear systems using MLF at switching instants:
 *  V_{σ(t_k)}(x(t_k)) ≤ V_{σ(t_{k-1})}(x(t_k)) */
bool dta_check_mlf_decrease(const DTA_SwitchedSystem* sys,
    const DTA_SwitchingSignal* sig, const double* x0,
    double (*V_func)(const double* x, int mode, int n));

#endif /* DTA_STABILITY_H */
