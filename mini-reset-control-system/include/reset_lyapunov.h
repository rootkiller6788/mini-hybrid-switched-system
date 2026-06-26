/* reset_lyapunov.h - Lyapunov-Based Stability Analysis for Reset Systems
 *
 * Provides Lyapunov function construction, passivity analysis, and
 * stability verification for reset control systems. Extends classical
 * Lyapunov theory to hybrid reset systems.
 *
 * Knowledge coverage:
 *   L2: Lyapunov stability for hybrid systems, passivity
 *   L4: Reset Lyapunov theorem, passivity theorem for reset systems,
 *       H_beta condition (Banos & Barreiro)
 *   L8: Quadratic Lyapunov functions for reset, dwell-time Lyapunov
 *
 * Ref: [BB12] Ch.5 on stability of reset systems;
 *      [NZT08] Nesic, Zaccarian, Teel (2008);
 *      Beker, Hollot, Chait, "Plant with integrator: an example of
 *      reset control overcoming limitations of linear feedback" (2004),
 *      IEEE TAC
 */

#ifndef RESET_LYAPUNOV_H
#define RESET_LYAPUNOV_H

#include "reset_core.h"
#include "reset_math.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L2: Lyapunov Function Types for Reset Systems
 * ================================================================ */

/** LyapFunction - quadratic Lyapunov function V(x) = x^T P x.
 *  P must be symmetric positive definite (P = P^T > 0). */
typedef struct {
    double    *P;       /**< Lyapunov matrix, n*n, symmetric */
    int        n;       /**< state dimension                 */
    bool       is_valid; /**< true if P is SPD               */
} LyapFunction;

/** LyapAnalysis - complete Lyapunov analysis result. */
typedef struct {
    LyapFunction  *V;            /**< the Lyapunov function found      */
    double         alpha;        /**< decay rate (V_dot <= -2*alpha*V) */
    double         beta;         /**< jump gain (V(x^+) <= beta*V(x))  */
    double         tau_dwell;    /**< required dwell time for stability */
    bool           is_stable;    /**< overall stability conclusion      */
    bool           is_quadratic; /**< true if quadratic Lyapunov works  */
    int            n_dim;        /**< system dimension                 */
    double         margin;       /**< stability margin (>0 if stable)  */
} LyapAnalysis;

/* ================================================================
 * L4: Fundamental Lyapunov Theorems for Reset Systems
 * ================================================================ */

/** Check positive definiteness of a symmetric matrix P.
 *  Uses Cholesky factorization attempt: P is SPD iff Cholesky succeeds.
 *  P: n*n symmetric matrix.
 *  Complexity: O(n^3).
 *  Returns true if P is symmetric positive definite. */
bool lyap_is_spd(int n, const double *P);

/** Construct a quadratic Lyapunov function for the base linear system
 *  by solving the Lyapunov equation:
 *
 *    A^T * P + P * A + Q = 0
 *
 *  where Q > 0 is chosen as identity. If A is Hurwitz, P > 0 exists
 *  and is unique.
 *
 *  Q: n*n positive definite matrix (if NULL, uses identity).
 *  P_out: n*n output Lyapunov matrix (pre-allocated).
 *
 *  Complexity: O(n^3) for Lyapunov equation.
 *  Returns true if a valid P > 0 was found. */
bool lyap_construct_linear(const ResetLinearBase *base,
                            const double *Q, double *P_out);

/** Verify the Lyapunov condition for a candidate P matrix:
 *    (1) P > 0 (positive definite)
 *    (2) A^T P + P A < 0 (negative definite, continuous strict decrease)
 *    (3) V(x^+) - V(x) <= 0 (non-increase at jumps)
 *
 *  This implements the hybrid Lyapunov theorem specialized to
 *  reset control systems.
 *
 *  Theorem (Banos & Barreiro, 2012, Theorem 5.1):
 *    If there exists P = P^T > 0 such that
 *      V_dot(x) = x^T (A^T P + P A) x < 0  for all x != 0 in C,
 *      Delta V = V(Ar x) - V(x) <= 0        for all x in D,
 *    then the reset system is globally asymptotically stable.
 *
 *  Returns true if all conditions are satisfied. */
bool lyap_verify_reset(const ResetLinearBase *base,
                        const ResetJumpMap *jump,
                        const double *P, double *margin);

/** Search for a Lyapunov function for a reset system.
 *
 *  Strategy:
 *    1. Solve Lyapunov equation for the base linear system: find P0.
 *    2. Check jump condition: P0 - Ar^T P0 Ar >= 0.
 *    3. If jump condition fails, try scaled P = alpha * P0.
 *    4. If still fails, solve LMI via iterative method.
 *
 *  The search is heuristic but effective for most reset systems.
 *
 *  Complexity: O(n^3 * n_trials).
 *  Returns a LyapAnalysis with results. Caller must free with
 *  lyap_analysis_free(). */
LyapAnalysis* lyap_search(const ResetSystem *rsys, int n_trials);

/** Free a LyapAnalysis result structure. */
void lyap_analysis_free(LyapAnalysis *analysis);

/* ================================================================
 * L4: Passivity Theorem for Reset Systems
 * ================================================================ */

/** Check if a linear system is strictly input passive:
 *
 *    V_dot <= u^T y - epsilon * ||u||^2  for some epsilon > 0
 *
 *  Equivalent to checking the KYP (Kalman-Yakubovich-Popov) lemma:
 *    There exists P = P^T > 0 such that
 *      [A^T P + P A,  P B - C^T;
 *       B^T P - C,   -(D + D^T) ] <= 0
 *
 *  Complexity: O(n^3). Returns true if strictly input passive. */
bool lyap_check_passivity(const ResetLinearBase *base, double epsilon);

/** Check if a reset system preserves passivity.
 *
 *  Theorem (Passivity of Reset Systems):
 *    If the base linear system is passive AND the reset map satisfies
 *      V(x^+) <= V(x^-)   (i.e., Ar^T P Ar <= P)
 *    where P is the storage function matrix from the KYP lemma,
 *    then the reset system is also passive.
 *
 *  This is important because passivity-based interconnections
 *  guarantee stability: feedback of passive systems is passive.
 *
 *  Ref: Carrasco, Banos, van der Schaft (2010),
 *       "A passivity-based approach to reset control systems stability"
 *
 *  Returns true if the reset system preserves passivity. */
bool lyap_check_reset_passivity(const ResetSystem *rsys);

/* ================================================================
 * L8: Advanced - Dwell-Time Lyapunov Conditions
 * ================================================================ */

/** Compute the minimum dwell time required for stability based
 *  on the Lyapunov function decay rate.
 *
 *  If V_dot <= -2*alpha*V during flow and V(x^+) <= beta*V(x) at jumps
 *  (with beta possibly > 1), then stability requires:
 *
 *    tau_dwell > ln(beta) / (2 * alpha)
 *
 *  This gives the minimum time between resets needed to guarantee
 *  that the overall hybrid trajectory is stable, even if individual
 *  jumps increase the Lyapunov function.
 *
 *  Ref: [NZT08] Theorem 2; Hespanha & Morse (1999) on dwell-time
 *
 *  Returns the required dwell time, or 0 if beta <= 1 (jumps always
 *  decrease V, so any dwell time works). */
double lyap_required_dwell_time(double alpha, double beta);

/** Evaluate the Lyapunov function V(x) = x^T P x at a given state.
 *  Complexity: O(n^2). */
double lyap_evaluate(const LyapFunction *V, const double *x);

/** Compute the directional derivative V_dot = gradV * f(x)
 *  for the flow vector field f(x) = A*x + B*u.
 *  Complexity: O(n^2). */
double lyap_flow_derivative(const LyapFunction *V,
                             const ResetLinearBase *base,
                             const double *x, const double *u);

/** Compute the jump increment Delta V = V(x^+) - V(x^-).
 *  x_after = Ar * x_before + Br * e.
 *  Complexity: O(n^2). */
double lyap_jump_increment(const LyapFunction *V,
                            const ResetJumpMap *jump,
                            const double *x_before, double e);

/* ================================================================
 * L4: Circle Criterion for Reset Systems
 * ================================================================ */

/** Verify the circle criterion for a reset control system.
 *
 *  For a reset system with FORE having reset ratio rho, the system
 *  is L2-stable if:
 *    (a) The base linear closed-loop is stable, and
 *    (b) || G_cl(jw) ||_inf * (1 + rho) / (1 - rho) < 1
 *        where G_cl is the complementary sensitivity.
 *
 *  This is a frequency-domain sufficient condition.
 *
 *  Ref: [BB12] Section 5.5; Zaccarian, Nesic, Teel (2005),
 *       "First order reset elements and the Clegg integrator revisited"
 *
 *  omega_min, omega_max: frequency range for Hinf search.
 *  n_points: number of frequency grid points.
 *  Returns true if circle criterion certifies stability. */
bool lyap_circle_criterion(const ResetLinearBase *plant,
                            const ResetSystem *controller,
                            double rho,
                            double omega_min, double omega_max,
                            int n_points);

/** Compute the H_beta matrix condition for reset stability.
 *
 *  H_beta(P) = [A^T P + P A + 2*beta*P,  P * Bp;
 *               Bp^T * P,                 0]
 *
 *  where (Ap, Bp, Cp) is the closed-loop system. If there exists
 *  P > 0 and beta > 0 such that H_beta(P) <= 0, then the reset
 *  system is asymptotically stable for any dwell time.
 *
 *  Returns the maximum beta for which the condition holds,
 *  or -1 if no such beta exists. */
double lyap_hbeta_condition(const ResetLinearBase *plant,
                             const ResetSystem *controller);

#ifdef __cplusplus
}
#endif

#endif /* RESET_LYAPUNOV_H */