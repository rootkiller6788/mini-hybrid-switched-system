#ifndef SWITCHED_STABILITY_H
#define SWITCHED_STABILITY_H

#include "switched_types.h"

/* ============================================================================
 * Switched Stability Analysis API (L1, L2, L4, L5)
 *
 * Implements fundamental stability analysis methods for switched systems:
 *   - Common Lyapunov Function (CLF) search
 *   - Multiple Lyapunov Functions (MLF) analysis
 *   - Lie-algebraic stability conditions
 *   - Dwell-time stability (Morse 1996)
 *   - Average dwell-time stability (Hespanha & Morse 1999)
 * ============================================================================ */

/* --------------------------------------------------------------------------
 * L1, L2: Stability Classification
 * -------------------------------------------------------------------------- */

/**
 * Classify stability of a switched system under a given switching signal.
 * Checks eigenvalues of each subsystem and applies the appropriate
 * stability criterion based on the switching signal type.
 *
 * Reference: Liberzon, D. (2003). Switching in Systems and Control. Ch. 2-3.
 * Complexity: O(p * n^3) for eigenvalue computation, p modes, dim n.
 */
SwitchedStabilityType sss_classify_stability(SwitchedSystem *sys);

/**
 * Check if the switched system is GUES under arbitrary switching using CLF.
 *
 * Condition: There exists P = P^T > 0 such that:
 *   A_i^T P + P A_i < 0  for all i in P
 *
 * If feasible, GUES follows with:
 *   ||x(t)|| <= sqrt(lambda_max(P)/lambda_min(P)) * exp(-gamma*t) * ||x(0)||
 *
 * Reference: Liberzon (2003), Theorem 2.1.
 */
bool sss_is_gues_arbitrary(SwitchedSystem *sys);

/**
 * Check stability under constrained switching with MLF.
 *
 * Uses the Branicky condition: for each mode i, V_i must decrease when
 * mode i is active, and the sequence {V_i(t_{i_k})} must be non-increasing
 * at switching instants.
 *
 * Reference: Branicky, M.S. (1998). IEEE TAC, 43(4), 475-482.
 */
bool sss_is_stable_mlf(SwitchedSystem *sys, MultipleLyapunovFunctions *mlf);

/* --------------------------------------------------------------------------
 * L4: Fundamental Theorems - Implementations
 * -------------------------------------------------------------------------- */

/**
 * Common Lyapunov Function Theorem (Liberzon 2003, Theorem 2.1).
 *
 * If all subsystems share a common Lyapunov function V(x) = x^T P x
 * with P = P^T > 0 and A_i^T P + P A_i < 0 for all i, then the switched
 * system is GUES under arbitrary switching.
 *
 * @param sys  Switched system with subsystems A_0,...,A_{p-1}
 * @param clf  Output: computed common Lyapunov function (if found)
 * @return     true if a CLF was found
 */
bool sss_clf_theorem(SwitchedSystem *sys, CommonLyapunovFunction *clf);

/**
 * Multiple Lyapunov Functions Theorem (Branicky 1998, Theorem 1).
 *
 * For each mode i, let V_i(x) = x^T P_i x be Lyapunov function when
 * mode i is active. The switched system is stable if:
 *   1. dV_i/dt < 0 when sigma(t) = i (active mode)
 *   2. V_{sigma(t_k)}(x(t_k)) <= V_{sigma(t_{k-1})}(x(t_k)) at switches
 *
 * Complexity: O(p * n^2) per check, where p = |P|.
 */
bool sss_mlf_theorem(SwitchedSystem *sys, MultipleLyapunovFunctions *mlf);

/**
 * Dwell-Time Stability Theorem (Morse 1996).
 *
 * If each subsystem A_i is Hurwitz (Re(lambda) <= -lambda_0 < 0),
 * then there exists tau_d > 0 such that dwell time >= tau_d => GUES.
 *
 * Minimum dwell time: tau_d = ln(K) / lambda_0,
 * where K depends on the Lyapunov function mismatch.
 *
 * Reference: Liberzon (2003), Section 3.2.1.
 */
void sss_dwell_time_theorem(SwitchedSystem *sys, DwellTimeAnalysis *dta);

/**
 * Average Dwell-Time Theorem (Hespanha & Morse 1999).
 *
 * Given lambda_0 > 0 and mu >= 1, if the switching signal has
 * average dwell time tau_a > ln(mu) / lambda_0, then GUES.
 *
 * N_sigma(T, t) <= N_0 + (T - t) / tau_a  for all T > t >= 0
 *
 * Reference: Liberzon (2003), Theorem 3.2.
 */
void sss_avg_dwell_theorem(SwitchedSystem *sys, DwellTimeAnalysis *dta);

/* --------------------------------------------------------------------------
 * L5: Algorithms - Lyapunov Function Computation
 * -------------------------------------------------------------------------- */

/**
 * Compute a Common Lyapunov Function (CLF) for a set of Hurwitz matrices.
 *
 * Uses iterative gradient descent on the LMI feasibility problem.
 * For commuting matrices (A_i A_j = A_j A_i), the CLF is explicit.
 *
 * Reference: Boyd et al. (1994). Linear Matrix Inequalities in System
 *   and Control Theory. SIAM. Section 2.2.
 *
 * Complexity: O(p * n^3 * max_iters), p = |P|, n = state_dim.
 */
bool sss_compute_clf(SwitchedSystem *sys, CommonLyapunovFunction *clf, int max_iters);

/**
 * Compute Multiple Lyapunov Functions (MLF) for each subsystem.
 *
 * For each Hurwitz A_i, solves Lyapunov equation:
 *   A_i^T P_i + P_i A_i = -Q  (default Q = I)
 *
 * Then checks MLF switching condition and computes mu parameter.
 */
void sss_compute_mlf(SwitchedSystem *sys, MultipleLyapunovFunctions *mlf);

/**
 * Solve continuous-time Lyapunov equation: A^T P + P A = -Q
 *
 * Uses vectorized form for n <= 4:
 *   (I kron A^T + A^T kron I) vec(P) = -vec(Q)
 *
 * Reference: Bartels & Stewart (1972). Comm. ACM, 15(9), 820-826.
 *
 * @param A  System matrix (n x n, must be Hurwitz)
 * @param Q  RHS matrix (n x n, positive definite)
 * @param P  Output solution matrix (n x n)
 * @param n  Dimension
 * @return   true if solution found
 */
bool sss_solve_lyapunov(const SwitchedMatrix *A, const SwitchedMatrix *Q,
                        SwitchedMatrix *P, int n);

/**
 * Check if a matrix is positive definite using Sylvester's criterion.
 * All leading principal minors must be > 0, and P must be symmetric.
 */
bool sss_is_positive_definite(const SwitchedMatrix *P);

/**
 * Compute mismatch parameter mu for MLF analysis.
 * mu = max_{i,j} (lambda_max(P_j) / lambda_min(P_i))
 * Appears in average dwell-time bound: tau_a* = ln(mu) / lambda_0.
 */
double sss_compute_mu(const MultipleLyapunovFunctions *mlf);

/**
 * Compute stability margin: lambda_0 = min_i -max_j Re(lambda_j(A_i)).
 * The slowest decay rate among all Hurwitz subsystems.
 */
double sss_compute_stability_margin(SwitchedSystem *sys);

/* --------------------------------------------------------------------------
 * L5: Algorithms - Lie-Algebraic Check
 * -------------------------------------------------------------------------- */

/**
 * Check the Lie-algebraic stability condition.
 *
 * If Lie{A_0, ..., A_{p-1}} is solvable, then the switched linear
 * system is GUES under arbitrary switching.
 *
 * Special cases:
 *   - All A_i commute pairwise -> nilpotent (strongest)
 *   - A_i are simultaneously triangularizable -> solvable
 *
 * Reference: Liberzon, Hespanha & Morse (1999). Systems & Control Letters.
 * Complexity: O(p^2 * n^3 * k) where k = derived length.
 */
bool sss_lie_algebraic_check(SwitchedSystem *sys, LieAlgebraCondition *la);

/* --------------------------------------------------------------------------
 * L5: Dwell-Time Computation
 * -------------------------------------------------------------------------- */

/**
 * Compute minimum dwell time: tau_d = ln(mu) / lambda_0.
 * If t_{k+1} - t_k >= tau_d for all k, the system is GUES.
 *
 * Reference: Liberzon (2003), Section 3.2.1.
 */
double sss_compute_min_dwell_time(double lambda_0, double mu);

/**
 * Compute required average dwell time: tau_a* = ln(mu) / lambda_0.
 *
 * Reference: Hespanha & Morse (1999), Liberzon (2003) Theorem 3.2.
 */
double sss_compute_avg_dwell_time(double lambda_0, double mu);

/**
 * Check if switch times satisfy dwell-time constraint.
 *
 * @param switch_times  Array of N+1 switch instants
 * @param n_switches    Number of switch intervals (= N)
 * @param tau_d         Minimum dwell time required
 * @return              true if all intervals >= tau_d
 */
bool sss_check_dwell_time(const double *switch_times, int n_switches, double tau_d);

/**
 * Check average dwell-time condition:
 *   N_sigma(T, 0) <= N_0 + T / tau_a for all T
 */
bool sss_check_avg_dwell(const SwitchingSignal *signal, double tau_a, int N0);

/**
 * Compute actual average dwell time: tau_avg = total_time / n_switches
 */
double sss_compute_actual_avg_dwell(const SwitchingSignal *signal);

/* --------------------------------------------------------------------------
 * L5: Simulation Methods
 * -------------------------------------------------------------------------- */

/**
 * Simulate switched system with Euler integration.
 * dx/dt = A_{sigma(t)} x(t)
 * Uses forward Euler: x(t+dt) = x(t) + dt * A_{sigma(t)} x(t)
 */
void sss_simulate_euler(SwitchedSystem *sys, const SolverConfig *cfg,
                        const SwitchSequence *seq);

/**
 * Simulate switched system with RK4 integration.
 * 4th-order Runge-Kutta for higher accuracy.
 *
 * Reference: Butcher (2008). Numerical Methods for ODEs. Wiley.
 */
void sss_simulate_rk4(SwitchedSystem *sys, const SolverConfig *cfg,
                      const SwitchSequence *seq);

/**
 * State-dependent switching simulation.
 * sigma(t) = phi(x(t)), where phi maps state to mode index.
 *
 * Example: sigma(x) = (x_1 > 0) ? 1 : 0 (bang-bang control)
 */
typedef int (*StateSwitchFunc)(const SwitchedVector *x, void *ctx);
void sss_simulate_state_dep(SwitchedSystem *sys, const SolverConfig *cfg,
                            StateSwitchFunc phi, void *phi_ctx);

#endif /* SWITCHED_STABILITY_H */
