#ifndef ETC_STABILITY_H
#define ETC_STABILITY_H

#include "etc_core.h"

/* ============================================================================
 * Event-Triggered Control — Stability Analysis
 *
 * Core stability concepts for event-triggered systems:
 *
 * 1. Input-to-State Stability (ISS) w.r.t. measurement error e(t)
 *    — The closed-loop is ISS from e to x if there exists a
 *      smooth Lyapunov function V(x) such that:
 *      V̇(x) ≤ −α(|x|) + γ(|e|)  with α, γ ∈ K∞.
 *
 * 2. Lyapunov-based triggering:
 *    — Γ(x,e) = γ(|e|) − σ α(|x|)  ⇒  V̇ ≤ −(1−σ)α(|x|)  when Γ < 0.
 *    — This enforces V̇ negative definite between events.
 *
 * 3. Zeno-freeness:
 *    — A strictly positive minimum inter-event time (MIET) exists iff
 *      the trigger function is Lipschitz in x and e.
 *
 * 4. L₂-gain bound:
 *    — The ETC system has finite L₂ gain from disturbance w to output z
 *      if a certain LMI condition holds.
 *
 * References:
 *   Khalil, H.K. (2002). Nonlinear Systems, 3rd ed. Prentice Hall.
 *   Sontag, E.D. (2008). Input to state stability. In Encyclopedia of
 *     Complexity and Systems Science. Springer.
 *   Tabuada, P. (2007). IEEE TAC 52(9).
 *   Heemels, W.P.M.H. et al. (2012). Intro to ETC & STC. IEEE TAC.
 * ============================================================================ */

/* --- Comparison Functions --- */

/** Evaluate a K-class function candidate at s.
 *  Common forms: α(s) = a·s (linear), α(s) = a·s² (quadratic),
 *  α(s) = a·s/(1+bs) (Michaelis-Menten saturation).
 *
 *  @param s    Non-negative input
 *  @param a    Scale parameter
 *  @param b    Shape parameter (used for saturating forms)
 *  @param form 0=linear, 1=quadratic, 2=saturating
 *  @return     α(s) */
double etc_kclass_function(double s, double a, double b, int form);

/** Verify that α is a valid K-class function on [0, s_max].
 *  Checks: α(0) = 0, α strictly increasing, α continuous.
 *  Samples at 100 points and verifies monotonicity.
 *
 *  @param alpha_fn  Function pointer
 *  @param a, b     Parameters
 *  @param form     Form index
 *  @param s_max    Upper bound of domain to check
 *  @return         True if all K-class properties hold on samples */
bool etc_verify_kclass(double (*alpha_fn)(double, double, double, int),
                        double a, double b, int form, double s_max);

/** Evaluate a KL-class function candidate β(s, t).
 *  Common form: β(s,t) = k · s · e^{−λ t}.
 *  Must satisfy: β(·,t) ∈ K for each fixed t, β(s,·) ↓ 0.
 *
 *  @param s     State norm (≥ 0)
 *  @param t     Time (≥ 0)
 *  @param k     Gain
 *  @param lam   Decay rate (λ > 0)
 *  @return      β(s, t) */
double etc_kl_function(double s, double t, double k, double lam);

/* --- ISS Lyapunov Analysis --- */

/** Check if V(x) = x^T P x is an ISS-Lyapunov function for the ETC system.
 *
 *  ISS condition: ∃ α₁, α₂ ∈ K∞, α₃, γ ∈ K such that
 *      α₁(|x|) ≤ V(x) ≤ α₂(|x|)
 *      V̇(x) ≤ −α₃(|x|) + γ(|e|)
 *
 *  For linear systems with quadratic V:
 *      α₁(s) = λ_min(P)·s²,  α₂(s) = λ_max(P)·s²
 *      V̇ = x^T (Acl^T P + P Acl) x + 2 x^T P B K e
 *         ≤ −λ_min(Q)|x|² + 2|PBK||x||e|
 *      α₃(s) = λ_min(Q)·s²,  γ(s) = (2|PBK|)²/(4θ)·s² + θ·s²
 *         = c₁ s²  with c₁ = |PBK|²/θ + θ  (for any θ > 0)
 *
 *  @param sys    ETC system with P and Acl properly set
 *  @param alpha1 Output: lower-bound K∞ function parameter (λ_min(P))
 *  @param alpha2 Output: upper-bound K∞ function parameter (λ_max(P))
 *  @param alpha3 Output: dissipation rate parameter
 *  @param gamma  Output: error-gain parameter
 *  @return       True if ISS-Lyapunov conditions hold */
bool etc_verify_iss_lyapunov(const ETCSystem* sys,
                              double* alpha1, double* alpha2,
                              double* alpha3, double* gamma);

/* --- Lyapunov Equation Solver --- */

/** Solve the continuous-time Lyapunov equation A^T P + P A = −Q.
 *  Uses Bartels-Stewart algorithm (n ≤ 16) or gradient descent (n > 16).
 *  Q is assumed to be identity (can be scaled).
 *
 *  @param A  Square matrix (n×n), must be Hurwitz
 *  @param P  Output: solution matrix (n×n), symmetric P > 0
 *  @param Q  Right-hand side (n×n), symmetric
 *  @param n  Dimension
 *  @return   True if solution found and P > 0 */
bool etc_solve_lyapunov(const ETCMatrix* A, ETCMatrix* P,
                         const ETCMatrix* Q, int n);

/** Compute the closed-loop Lyapunov function for (A+BK).
 *  Solves (A+BK)^T P + P (A+BK) = −I and stores result in sys->V.
 *
 *  @param sys  ETC system (A, B, K must be set)
 *  @return     True if P > 0 found (closed-loop is Hurwitz) */
bool etc_compute_clf(ETCSystem* sys);

/* --- L2-Gain Analysis --- */

/** Compute L₂-gain bound from disturbance w to state x for ETC.
 *
 *  ẋ = Acl x + Bw w  (Bw: disturbance input matrix)
 *  Under event-triggering with threshold σ:
 *     ||x||_{L₂[0,T]} ≤ γ ||w||_{L₂[0,T]} + β(|x(0)|, T)
 *
 *  The L₂-gain γ is bounded by solving the LMI:
 *     [Acl^T P + P Acl + I + (1/σ) P BK K^T B^T P    P Bw]
 *     [Bw^T P                                       −γ² I]  < 0
 *
 *  @param sys   ETC system with Acl, B, K set
 *  @param Bw    Disturbance input matrix (n × nw)
 *  @param sigma Trigger threshold
 *  @param gamma Output: L₂-gain upper bound
 *  @return      True if finite gain found */
bool etc_compute_l2_gain(const ETCSystem* sys, const ETCMatrix* Bw,
                          double sigma, double* gamma);

/* --- Zeno Analysis --- */

/** Determine if the ETC system is Zeno-free.
 *
 *  Zeno behavior: lim_{k→∞} t_k = T_∞ < ∞ (infinite events in finite time).
 *
 *  For static triggering with Γ = |e| − σ|x|, Zeno-freeness is guaranteed
 *  if the system dynamics satisfy a Lipschitz condition:
 *     d/dt (|e|/|x|) ≤ L (1 + |e|/|x|) for some L < ∞.
 *
 *  This function checks the theoretical Zeno-free condition and also
 *  detects Zeno from simulation data.
 *
 *  @param sys  ETC system (with event history populated)
 *  @return     True if Zeno-free (provably or empirically) */
bool etc_check_zeno_free(const ETCSystem* sys);

/** Compute the theoretical Zeno time if system is Zeno.
 *  For a linear system with static trigger, if σ > σ_max, V̇ may be
 *  positive, potentially leading to Zeno.
 *
 *  @param sys  ETC system
 *  @return     Estimated Zeno time, or INFINITY if Zeno-free */
double etc_compute_zeno_time(const ETCSystem* sys);

/* --- Practical Stability --- */

/** Verify practical stability of the ETC system.
 *
 *  A system is practically stable if trajectories converge to a
 *  bounded set B(r) = {x : |x| ≤ r} rather than the origin.
 *  This occurs when σ is too large or when using absolute thresholds.
 *
 *  @param sys        ETC system (after simulation)
 *  @param radius     Output: estimated convergence radius r
 *  @return           True if practically stable */
bool etc_check_practical_stability(const ETCSystem* sys, double* radius);

/** Compute ultimate bound for state under static event-triggering.
 *  Ultimate bound: |x(t)| ≤ b eventually, where
 *     b = σ · ||BK|| · λ_max(P) / ( (1−σ) λ_min(Q) ) · sup_t |x(t)|.
 *
 *  @param sys  ETC system with P, Acl, BK set
 *  @param sigma Trigger threshold
 *  @return     Ultimate bound b */
double etc_ultimate_bound(const ETCSystem* sys, double sigma);

#endif /* ETC_STABILITY_H */
