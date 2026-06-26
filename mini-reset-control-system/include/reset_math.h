/* reset_math.h - Mathematical Utilities for Reset Control Analysis
 *
 * Provides numerical linear algebra primitives needed for reset
 * control system analysis: matrix operations, eigenvalue computation,
 * Lyapunov equation solver, Riccati equation, and transfer function
 * evaluation at complex frequencies.
 *
 * Knowledge coverage:
 *   L3: Matrix algebra, eigenvalue decomposition
 *   L4: Lyapunov stability equation (AX + XA^T + Q = 0)
 *   L5: QR algorithm, Bartels-Stewart algorithm, matrix sign function
 *
 * Ref: Golub & Van Loan, "Matrix Computations" (2013), 4th ed.
 *      Bartels & Stewart, "Algorithm 432: Solution of the matrix
 *      equation AX + XB = C" (1972), Comm. ACM
 */

#ifndef RESET_MATH_H
#define RESET_MATH_H

#include "reset_core.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L3: Basic Matrix Operations
 * ================================================================ */

/** Matrix-Vector Multiply: y = alpha * A*x + beta * y.
 *  A in R^{m*n}, x in R^n, y in R^m.
 *  alpha, beta: scalars (use 1,0 for pure multiply).
 *  Complexity: O(m*n). */
void mat_vec_mul(int m, int n, const double *A, const double *x,
                  double *y, double alpha, double beta);

/** Matrix-Matrix Multiply: C = alpha * A*B + beta * C.
 *  A in R^{m*k}, B in R^{k*n}, C in R^{m*n}.
 *  All matrices stored row-major.
 *  Complexity: O(m*k*n). */
void mat_mat_mul(int m, int k, int n,
                  const double *A, const double *B, double *C,
                  double alpha, double beta);

/** Matrix transpose in-place: B = A^T.
 *  A is m*n, B is n*m (pre-allocated).
 *  Complexity: O(m*n). */
void mat_transpose(int m, int n, const double *A, double *B);

/** Matrix inversion via LU decomposition with partial pivoting.
 *  A: input matrix (n*n), overwritten with LU factors.
 *  A_inv: output inverse matrix (n*n), pre-allocated.
 *  Complexity: O(n^3). Returns 0 on success, -1 if singular. */
int mat_invert(int n, const double *A, double *A_inv);

/** Solve linear system A * x = b via LU decomposition.
 *  A: n*n matrix (not overwritten). b: n-vector (overwritten with x).
 *  Complexity: O(n^3). Returns 0 on success, -1 if singular. */
int mat_solve(int n, const double *A, double *bx);

/** Compute the determinant of an n*n matrix via LU decomposition.
 *  Complexity: O(n^3). Returns 0.0 if singular (within eps). */
double mat_determinant(int n, const double *A);

/** Compute the trace (sum of diagonal elements) of an n*n matrix.
 *  Complexity: O(n). */
double mat_trace(int n, const double *A);

/** Compute the Frobenius norm: sqrt(sum_{i,j} A_{ij}^2).
 *  Complexity: O(m*n). */
double mat_frobenius_norm(int m, int n, const double *A);

/** Identity matrix: sets A (n*n) to identity.
 *  Complexity: O(n^2). */
void mat_set_identity(int n, double *A);

/** Scale matrix in place: A = alpha * A.
 *  Complexity: O(m*n). */
void mat_scale(int m, int n, double *A, double alpha);

/** Copy matrix: dst = src.
 *  Complexity: O(m*n). */
void mat_copy(int m, int n, const double *src, double *dst);

/** Zero matrix: sets all entries to 0.
 *  Complexity: O(m*n). */
void mat_zero(int m, int n, double *A);

/* ================================================================
 * L5: Eigenvalue Computation - QR Algorithm
 * ================================================================ */

/** Compute eigenvalues of a real square matrix via the
 *  Francis double-shift QR algorithm with Hessenberg reduction.
 *
 *  A: n*n matrix (overwritten during computation).
 *  eigen_real, eigen_imag: pre-allocated arrays of length n.
 *  max_iter: maximum QR iterations (typically 100*n).
 *  tol: convergence tolerance.
 *
 *  Complexity: O(n^3) average case.
 *  Returns the number of converged eigenvalues. */
int mat_eigenvalues_qr(int n, double *A, double *eigen_real,
                        double *eigen_imag, int max_iter, double tol);

/** Compute eigenvalues using the power iteration method
 *  (only dominant eigenvalue). Good for large sparse matrices.
 *  A: n*n matrix. x: initial guess and output eigenvector (length n).
 *  lambda: output dominant eigenvalue.
 *  max_iter: max iterations. tol: convergence tolerance.
 *  Complexity: O(n^2 * iter). Returns number of iterations. */
int mat_power_iteration(int n, const double *A, double *x,
                         double *lambda, int max_iter, double tol);

/* ================================================================
 * L4: Lyapunov Equation Solver - Bartels-Stewart Algorithm
 * ================================================================ */

/** Solve the continuous-time Lyapunov equation:
 *
 *    A * X + X * A^T + Q = 0
 *
 *  where A in R^{n*n}, Q = Q^T in R^{n*n}, X = X^T in R^{n*n}.
 *
 *  Uses the Bartels-Stewart algorithm [BS72]:
 *    1. Reduce A to real Schur form: T = U^T * A * U
 *    2. Transform Q: Q_tilde = U^T * Q * U
 *    3. Solve T*Y + Y*T^T + Q_tilde = 0 (back-substitution)
 *    4. Recover X = U * Y * U^T
 *
 *  Prerequisite: A must be stable (all eigenvalues in left half-plane)
 *  for a unique solution to exist.
 *
 *  Complexity: O(n^3).
 *  Returns 0 on success, -1 if A is not stable or singular. */
int mat_lyapunov_solve(int n, const double *A, const double *Q, double *X);

/** Solve the discrete-time Lyapunov equation:
 *
 *    A * X * A^T - X + Q = 0
 *
 *  Complexity: O(n^3). Returns 0 on success. */
int mat_dlyap_solve(int n, const double *A, const double *Q, double *X);

/* ================================================================
 * L5: Algebraic Riccati Equation (ARE)
 * ================================================================ */

/** Solve the continuous-time algebraic Riccati equation:
 *
 *    A^T * P + P * A - P * B * R^{-1} * B^T * P + Q = 0
 *
 *  via the Hamiltonian matrix method.
 *
 *  H = [A,  -B*R^{-1}*B^T;
 *       -Q, -A^T          ]
 *
 *  The stable invariant subspace of H gives the solution P.
 *
 *  Complexity: O(n^3). Returns 0 on success, -1 if no stabilizing solution.
 *
 *  Ref: Laub (1979), "A Schur method for solving algebraic Riccati equations"
 */
int mat_riccati_solve(int n, int m,
                       const double *A, const double *B,
                       const double *Q, const double *R,
                       double *P);

/* ================================================================
 * L3: Transfer Function Evaluation
 * ================================================================ */

/** Evaluate the transfer function G(s) = C*(sI - A)^{-1}*B + D
 *  at a complex frequency s = sigma + j*omega.
 *
 *  G_real, G_imag: output real and imaginary parts of the transfer
 *  function matrix (size p*m, row-major).
 *
 *  Complexity: O(n^3) for matrix inversion.
 *  Returns 0 on success, -1 on near-singular (sI-A). */
int mat_transfer_func(int n, int m, int p,
                       const double *A, const double *B,
                       const double *C, const double *D,
                       double sigma, double omega,
                       double *G_real, double *G_imag);

/** Evaluate the frequency response magnitude at s = j*omega.
 *  ||G(jw)||_2 for each input-output pair.
 *  mag: output array of size p*m.
 *  Complexity: O(n^3). */
int mat_freqresp_mag(int n, int m, int p,
                      const double *A, const double *B,
                      const double *C, const double *D,
                      double omega, double *mag);

/** Compute the maximum singular value sigma_max(G(jw)).
 *  Uses power iteration on G(jw)^H * G(jw).
 *  Complexity: O(n^3 + p*m*iter). */
double mat_sigma_max(int n, int m, int p,
                      const double *A, const double *B,
                      const double *C, const double *D,
                      double omega, int max_iter, double tol);

#ifdef __cplusplus
}
#endif

#endif /* RESET_MATH_H */