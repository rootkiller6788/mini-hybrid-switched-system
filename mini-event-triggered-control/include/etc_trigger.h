#ifndef ETC_TRIGGER_H
#define ETC_TRIGGER_H

#include "etc_core.h"

/* ============================================================================
 * Event-Triggered Control — Trigger Conditions & Threshold Functions
 *
 * The triggering condition determines when a new control update is needed.
 * Three canonical families are implemented:
 *
 * 1. Static (state-relative):   Γ(x,e) = |e| − σ|x|          [Tabuada 2007]
 * 2. Quadratic (energy-based):  Γ(x,e) = |e|² − σ|x|²       [Heemels 2012]
 * 3. Mixed (time-dependent):    Γ(x,e,t) = |e|² − σ|x|² − ε e^{-αt}  [Wang 2015]
 * 4. Absolute:                  Γ(x,e) = |e| − ε            [Astrom 2002]
 * 5. Dynamic:            Γ(x,e,η) = θ η + σ|x|² − |e|²     [Girard 2015]
 *
 * Reference:
 *   Tabuada, P. (2007). Event-triggered real-time scheduling of stabilizing
 *     control tasks. IEEE TAC, 52(9), 1680–1685.
 *   Girard, A. (2015). Dynamic triggering mechanisms for event-triggered control.
 *     IEEE TAC, 60(7), 1992–1997.
 *   Wang, X. & Lemmon, M.D. (2011). Event-triggering in distributed NCS.
 *     IEEE TAC, 56(3), 586–601.
 * ============================================================================ */

/* --- Threshold Parameter Structures --- */

/** Static threshold parameters.
 *  Event when |e| > σ|x|. Smaller σ ⇒ more events, better performance. */
typedef struct {
    double sigma;            /* Proportional threshold (0 < σ < 1) */
    double sigma_min;        /* Minimum admissible σ (from stability condition) */
    bool is_admissible;      /* True if σ satisfies stability condition */
} ETCStaticThreshold;

/** Quadratic threshold parameters.
 *  Event when |e|² > σ|x|². Equivalent to static with squared norms. */
typedef struct {
    double sigma;            /* Quadratic threshold (0 < σ < 1) */
    double sigma_critical;   /* σ* = λ_min(Q) / (2||PBK||)  (critical value) */
    bool is_stabilizing;     /* True if σ < σ_critical */
} ETCQuadraticThreshold;

/** Mixed threshold parameters (time-varying bound).
 *  Event when |e|² > σ|x|² + ε e^{-α t}.
 *  The decaying exponential term provides robustness to initial transients. */
typedef struct {
    double sigma;            /* Proportional part */
    double epsilon;          /* Magnitude of exponential decay term */
    double alpha;            /* Decay rate (α > 0) */
    double t0;               /* Reference time for decay */
} ETCMixedThreshold;

/** Absolute threshold parameters (constant bound).
 *  Event when |e| > ε. Simpler but no stability guarantees in general. */
typedef struct {
    double epsilon;          /* Absolute error tolerance */
    double epsilon_min;      /* Minimum epsilon for practical stability */
    bool is_practically_stable;
} ETCAbsoluteThreshold;

/** Dynamic threshold with internal variable η.
 *  η̇ = −β η + σ|x|² − |e|², η(0) = η₀ ≥ 0.
 *  Event when η(t) + θ(σ|x|² − |e|²) ≤ 0.
 *  Guarantees: η(t) ≥ 0 for all t, fewer events than static counterpart. */
typedef struct {
    double theta;            /* Dynamic weight (0 < θ < 1) */
    double beta;             /* Decay rate for η (β > 0) */
    double eta0;             /* Initial η value (≥ 0) */
    double sigma;            /* Underlying static threshold */
    bool is_monotone;        /* True if η(t) is monotone non-increasing */
} ETCDynamicThreshold;

/** Self-triggered threshold parameters.
 *  Next event time t_{k+1} = t_k + τ(x(t_k)) computed at event k.
 *  τ(x) = min{ t > 0 : |e^{Acl t} x − x| = σ|x| } under zero-order hold. */
typedef struct {
    double sigma;            /* Relative threshold used to compute τ */
    double tau_min;          /* Minimum allowed τ (MIET guarantee) */
    double tau_max;          /* Maximum allowed τ */
    double prediction_horizon;
} ETCSelfTriggeredThreshold;

/* --- Trigger Function Library --- */

/** Standard static trigger: Γ(x,e) = |e| − σ|x|.
 *  Returns positive value when event should fire. */
double etc_trigger_static(const ETCVector* x, const ETCVector* e,
                           double sigma, double epsilon);

/** Quadratic trigger: Γ(x,e) = |e|² − σ|x|².
 *  Equivalent to static but with squared norms — smoother threshold surface. */
double etc_trigger_quadratic(const ETCVector* x, const ETCVector* e,
                              double sigma, double epsilon);

/** Mixed trigger: Γ(x,e) = |e|² − σ|x|² − ε e^{-α(t−t0)}.
 *  Note: t, t0 are passed implicitly through epsilon and sigma adjustments;
 *  the caller must update epsilon_effective = ε * exp(−α*(t−t0)). */
double etc_trigger_mixed(const ETCVector* x, const ETCVector* e,
                          double sigma, double epsilon_effective);

/** Absolute trigger: Γ(x,e) = |e| − ε.
 *  Simplest form; no state-dependence guarantees Zeno-freeness. */
double etc_trigger_absolute(const ETCVector* x, const ETCVector* e,
                             double sigma, double epsilon);

/** Dynamic trigger: Γ(x,e,η) = θ η + σ|x|² − |e|².
 *  Requires eta to be maintained externally (dynamic variable). */
double etc_trigger_dynamic(const ETCVector* x, const ETCVector* e,
                            double eta, double sigma, double theta);

/** Self-triggered next-interval computation.
 *  Given current sampled state x_k and closed-loop matrix Acl,
 *  estimates the time until |e(t)| = σ|x(t)| under zero-order hold u ≡ K x_k.
 *  Uses: x(t) = e^{Acl t} x_k + ∫₀ᵗ e^{Acl(t−s)} B K x_k ds.
 *
 *  @param x_k        Sampled state at event
 *  @param Acl        Closed-loop matrix (A+BK), n×n
 *  @param BK         Precomputed BK product matrix, n×n
 *  @param sigma      Threshold
 *  @param tau_min    Hard lower bound on returned time
 *  @param tau_max    Hard upper bound on returned time
 *  @return           Estimated next inter-event time τ */
double etc_self_triggered_interval(const ETCVector* x_k,
                                    const ETCMatrix* Acl,
                                    const ETCMatrix* BK,
                                    double sigma,
                                    double tau_min,
                                    double tau_max);

/* --- Trigger Analysis --- */

/** Compute the maximum admissible σ for static triggering.
 *  From stability condition: σ_max = λ_min(Q) / (2||PBK||).
 *  Requires Lyapunov matrix P and closed-loop matrices.
 *
 *  @param P   Lyapunov matrix (from Acl^T P + P Acl = −Q)
 *  @param BK  Product matrix B × K
 *  @return    σ_max — any σ ∈ (0, σ_max) guarantees asymptotic stability */
double etc_compute_sigma_max(const ETCMatrix* P, const ETCMatrix* BK);

/** Check if a given σ is stabilizing for static ETC.
 *  Verifies that σ < λ_min(Q) / (2||PBK||). */
bool etc_is_sigma_stabilizing(const ETCMatrix* P, const ETCMatrix* BK,
                               double sigma);

/** Compute the minimum inter-event time lower bound for static ETC.
 *  τ_min = σ / (||Acl|| + σ ||BK||)  [Tabuada 2007, Lemma 4].
 *
 *  @param Acl   Closed-loop matrix
 *  @param BK    Product B×K
 *  @param sigma Threshold
 *  @return      Positive lower bound on inter-event times */
double etc_compute_iet_lower_bound(const ETCMatrix* Acl, const ETCMatrix* BK,
                                    double sigma);

/** Compute the dynamic trigger ODE right-hand side.
 *  η̇ = −β η + σ|x|² − |e|².
 *
 *  @param eta   Current η
 *  @param beta  Decay rate
 *  @param x     Current state
 *  @param e     Measurement error
 *  @param sigma Threshold
 *  @return      dη/dt */
double etc_dynamic_trigger_rhs(double eta, double beta,
                                const ETCVector* x, const ETCVector* e,
                                double sigma);

/** Design a static threshold pair (σ, ε) for a given performance target.
 *  Balances inter-event time (communication reduction) against
 *  performance degradation (larger error bound).
 *
 *  @param desired_iet  Target average inter-event time
 *  @param performance  1.0 = no degradation, 0.0 = loose triggering
 *  @param Acl          Closed-loop matrix for bound computation
 *  @param BK           Product B×K
 *  @param sigma_out    Output: designed σ
 *  @param epsilon_out  Output: designed ε (for mixed trigger) */
void etc_design_threshold(double desired_iet, double performance,
                           const ETCMatrix* Acl, const ETCMatrix* BK,
                           double* sigma_out, double* epsilon_out);

#endif /* ETC_TRIGGER_H */
