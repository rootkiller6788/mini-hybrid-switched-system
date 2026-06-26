#ifndef SWITCHED_LYAPUNOV_H
#define SWITCHED_LYAPUNOV_H

#include "switched_types.h"

/* ============================================================================
 * Lyapunov Function Computation for Switched Systems (L3, L4, L5)
 *
 * Implements:
 *   - Common Lyapunov function matrix P via LMI methods
 *   - Multiple Lyapunov functions via Lyapunov equation solving
 *   - Eigenvalue bounds for P matrices
 *   - Lyapunov function evaluation and derivative computation
 * ============================================================================ */

/* --------------------------------------------------------------------------
 * L3: Lyapunov Matrix Computation
 * -------------------------------------------------------------------------- */

/**
 * Compute eigenvalues of a symmetric 2x2 matrix.
 * For P = [[a, b], [b, c]]:
 *   lambda = (a+c)/2 +/- sqrt(((a-c)/2)^2 + b^2)
 *
 * Complexity: O(1), closed-form solution.
 */
void sss_eig_sym_2x2(const SwitchedMatrix *P, double *eig1, double *eig2);

/**
 * Compute min and max eigenvalues of a small symmetric matrix (n <= 4).
 * Uses power iteration for largest, inverse power for smallest.
 *
 * Reference: Golub & Van Loan (2013). Matrix Computations. 4th ed.
 * Complexity: O(n^2 * iters) per eigenvalue.
 */
void sss_eig_sym(const SwitchedMatrix *P, double *eig_min, double *eig_max, int max_iters);

/**
 * Compute eigenvalues of a 2x2 real matrix (not necessarily symmetric).
 *
 * For A = [[a, b], [c, d]]:
 *   lambda = (tr(A) +/- sqrt(tr(A)^2 - 4*det(A))) / 2
 */
void sss_eig_2x2(const SwitchedMatrix *A, EigenvalueResult *e1, EigenvalueResult *e2);

/**
 * Compute spectral radius rho(A) = max_i |lambda_i(A)|.
 * Uses QR algorithm for eigenvalue computation.
 *
 * Complexity: O(n^3 * iters) for full eigenvalue decomposition.
 */
double sss_spectral_radius(const SwitchedMatrix *A);

/* --------------------------------------------------------------------------
 * L4: Lyapunov Equation Solvers
 * -------------------------------------------------------------------------- */

/**
 * Solve Lyapunov equation A^T P + P A = -Q for n=2 using direct formula.
 *
 * For A = [[a, b], [c, d]], Q = I, the solution is:
 *   P = (1/(2*det(A)*tr(A))) * [[-det(A)-c^2-d^2, ab+cd],
 *                                 [ab+cd, -det(A)-a^2-b^2]]
 *
 * Requires A to be Hurwitz (tr(A) < 0, det(A) > 0).
 *
 * Reference: Khalil (2002). Nonlinear Systems. 3rd ed., Section 3.3.
 * Complexity: O(1).
 */
bool sss_lyap_2x2(const SwitchedMatrix *A, SwitchedMatrix *P);

/**
 * Solve Lyapunov equation A^T P + P A = -Q for n=3 using Kronecker method.
 *
 * Constructs the 9x9 Kronecker system:
 *   (I kron A^T + A^T kron I) vec(P) = -vec(Q)
 * and solves via Gaussian elimination.
 *
 * Complexity: O(n^6) for Kronecker construction, O(n^6) for solve.
 */
bool sss_lyap_3x3(const SwitchedMatrix *A, const SwitchedMatrix *Q, SwitchedMatrix *P);

/**
 * General Lyapunov equation solver for n <= 6 using vectorization.
 *
 * Uses Kronecker product method: (I kron A^T + A^T kron I) vec(P) = -vec(Q).
 * Solves the resulting n^2 x n^2 linear system via Gaussian elimination
 * with partial pivoting.
 *
 * Reference: Bartels & Stewart (1972). Comm. ACM, 15(9), 820-826.
 * Complexity: O(n^6).
 *
 * @param A  System matrix (Hurwitz required)
 * @param Q  Right-hand side (positive definite by default)
 * @param P  Output Lyapunov matrix
 * @param n  Dimension (n <= 6)
 * @return   true if successful (A is Hurwitz, solution exists)
 */
bool sss_lyap_solve(const SwitchedMatrix *A, const SwitchedMatrix *Q, SwitchedMatrix *P, int n);

/* --------------------------------------------------------------------------
 * L5: Common Lyapunov Function Search
 * -------------------------------------------------------------------------- */

/**
 * Check if a given P satisfies the CLF condition for all subsystems.
 *
 * Condition: P = P^T > 0 and A_i^T P + P A_i < 0 for all i.
 *
 * This is the core verification step for the CLF theorem.
 */
bool sss_clf_verify(const SwitchedMatrix *P, SwitchedSubsystem **modes, int n_modes);

/**
 * Compute CLF via projected gradient descent on the LMI cone.
 *
 * Minimizes: f(P) = max_i lambda_max(A_i^T P + P A_i)
 * subject to: P = P^T, lambda_min(P) >= eps > 0
 *
 * Algorithm: iterate P_{k+1} = P_k - alpha * grad_f(P_k),
 * then project back onto the positive definite cone.
 *
 * Reference: Boyd et al. (1994). LMI in System and Control Theory. SIAM.
 * Complexity: O(max_iters * p * n^3).
 *
 * @param modes      Array of subsystems
 * @param n_modes    Number of modes
 * @param n          State dimension
 * @param clf        Output CLF structure
 * @param max_iters  Maximum gradient descent iterations
 * @return           true if a feasible P was found
 */
bool sss_clf_gradient_descent(SwitchedSubsystem **modes, int n_modes, int n,
                               CommonLyapunovFunction *clf, int max_iters);

/* --------------------------------------------------------------------------
 * L5: Multiple Lyapunov Function Computation
 * -------------------------------------------------------------------------- */

/**
 * Compute MLF by solving Lyapunov equation for each Hurwitz subsystem.
 *
 * For each mode i with Hurwitz A_i, solves:
 *   A_i^T P_i + P_i A_i = -I
 *
 * Then computes the switching attenuation parameter mu.
 */
void sss_mlf_compute_all(SwitchedSystem *sys, MultipleLyapunovFunctions *mlf);

/**
 * Verify the MLF non-increasing sequence condition at switch instants.
 *
 * For each switch from mode i to mode j at time t_k:
 *   V_j(x(t_k)) <= mu * V_i(x(t_k))
 * must hold for some mu >= 1.
 *
 * This implements the Branicky condition verification.
 */
bool sss_mlf_verify_sequence(const SwitchedSystem *sys, const MultipleLyapunovFunctions *mlf,
                              const SwitchingSignal *signal);

/**
 * Estimate Lyapunov function value: V(x) = x^T P x.
 *
 * @param P  Lyapunov matrix
 * @param x  State vector
 * @return   V(x) = x^T P x
 */
double sss_lyap_eval(const SwitchedMatrix *P, const SwitchedVector *x);

/**
 * Estimate Lyapunov derivative: dV/dt = x^T (A^T P + P A) x.
 *
 * For mode i with dynamics dx/dt = A_i x:
 *   dV_i/dt(x) = x^T (A_i^T P_i + P_i A_i) x
 *
 * If A_i^T P_i + P_i A_i = -Q_i < 0, then dV_i/dt < 0 for x != 0.
 */
double sss_lyap_derivative(const SwitchedMatrix *P, const SwitchedMatrix *A,
                            const SwitchedVector *x);

#endif /* SWITCHED_LYAPUNOV_H */
