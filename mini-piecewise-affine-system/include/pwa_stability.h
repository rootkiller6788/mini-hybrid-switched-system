/**
 * @file pwa_stability.h
 * @brief Stability Analysis for PWA Systems — L4 Fundamental Theorems
 *
 * Stability analysis for piecewise affine systems uses piecewise
 * quadratic (PWQ) Lyapunov functions and the S-procedure to
 * handle the polyhedral partition.
 *
 * A PWQ Lyapunov function has the form:
 *   V(x) = x^T P_i x + 2 q_i^T x + r_i   for x ∈ R_i
 *
 * where P_i = P_i^T is a symmetric matrix for each region.
 *
 * Stability conditions (CT-PWA):
 *   V(x) > 0  for all x ≠ 0 in each R_i (positivity)
 *   dV/dt(x) < 0 for all x ≠ 0 in each R_i (decrease)
 *   V is continuous across region boundaries
 *
 * The S-procedure converts these conditions to Linear Matrix
 * Inequalities (LMIs) that can be solved via semidefinite programming.
 *
 * References:
 *   Johansson, M. & Rantzer, A. (1998). "Computation of piecewise
 *     quadratic Lyapunov functions for hybrid systems."
 *     IEEE TAC, 43(4):555-559.
 *   Boyd, S., El Ghaoui, L., Feron, E., & Balakrishnan, V. (1994).
 *     "Linear Matrix Inequalities in System and Control Theory."
 *     SIAM Studies in Applied Mathematics, Vol. 15.
 *   Pettersson, S. & Lennartson, B. (1999). "Exponential stability
 *     of hybrid systems using piecewise quadratic Lyapunov functions
 *     resulting in LMIs." IFAC World Congress.
 *
 * Knowledge coverage:
 *   L4: PWQ Lyapunov stability, S-procedure, invariant sets,
 *       LMI formulation for PWA stability
 *
 * Nine-school course alignment:
 *   MIT 6.241J — Lec 13 Lyapunov stability
 *   Stanford AA203 — Lec 6 Lyapunov theory
 *   Berkeley EE222 — Ch 4 Lyapunov stability
 *   CMU 24-677 — Lec 15 Switched system stability
 *   Princeton MAE 546 — Lec 12 Lyapunov methods
 *   Caltech CDS140 — Lec 8 Stability theory
 *   Cambridge 4F3 — Lec 4 Lyapunov for nonlinear
 *   Oxford B4 — Lec 7 Stability analysis
 *   ETH 227-0216 — Lec 14 Stability
 */

#ifndef PWA_STABILITY_H
#define PWA_STABILITY_H

#include "pwa_defs.h"
#include "pwa_geometry.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
 * L4: Piecewise Quadratic Lyapunov Functions
 *===========================================================================*/

/**
 * @brief Piecewise Quadratic (PWQ) Lyapunov function.
 *
 * V(x) = x^T P_i x + 2 q_i^T x + r_i  for x ∈ R_i
 *
 * Each region has its own quadratic form.
 * Continuity condition: V_i(x) = V_j(x) for x on boundary R_i ∩ R_j.
 */
typedef struct {
    int      n_regions;     /**< Number of regions */
    int      n_state;       /**< State dimension */
    double **P;             /**< P_i matrices (n_state × n_state, symmetric) */
    double **q;             /**< q_i vectors (n_state) */
    double  *r;             /**< r_i scalars */
    int      is_continuous; /**< 1 if V is continuous across boundaries */
    int      is_positive;   /**< 1 if V(x) > 0 for all x ≠ 0 */
    int      is_decreasing; /**< 1 if dV/dt < 0 along trajectories */
} PWQLyapunov;

/**
 * @brief S-procedure parameters for PWA stability analysis.
 *
 * The S-procedure relaxes a condition "f_0(x) ≥ 0 whenever
 * f_i(x) ≥ 0, i=1..m" to: ∃ τ_i ≥ 0 s.t. f_0(x) - Σ τ_i f_i(x) ≥ 0 ∀x.
 *
 * For PWA systems, f_i are the region constraints H_i x ≤ K_i.
 */
typedef struct {
    int      n_state;       /**< State dimension */
    int      n_regions;     /**< Number of regions */
    double   eps;           /**< Tolerance for strict inequalities */
    double **T;             /**< S-procedure multipliers for each region */
    double  *alpha;         /**< Decay rate α: dV/dt ≤ -α V */
} PWASProcedure;

/**
 * @brief Linear Matrix Inequality (LMI) for PWA stability.
 *
 * An LMI has the form: F_0 + Σ x_i F_i ≽ 0 (positive semidefinite)
 *
 * For PWA stability, the LMI encodes the PWQ Lyapunov conditions.
 */
typedef struct {
    int      n_vars;        /**< Number of decision variables */
    int      lmi_dim;       /**< Dimension of the LMI matrices */
    int      n_lmis;        /**< Number of LMI constraints */
    double **F;             /**< LMI matrices F_0, F_1, ..., F_{n_vars} */
} PWALMI;

/**
 * @brief Invariant set for a PWA system.
 *
 * A set S is invariant if x(0) ∈ S ⇒ x(t) ∈ S for all t ≥ 0.
 * A set S is positively invariant if x ∈ S ⇒ f(x) ∈ S (DT) or
 * the vector field points inward (CT).
 */
typedef struct {
    int             n_state;     /**< State dimension */
    PWAPolyhedron   poly;        /**< Polyhedral invariant set */
    int             region_id;   /**< Region containing this set */
    int             is_invariant; /**< Verification flag */
} PWAInvariantSet;

/*===========================================================================
 * L4: PWQ Lyapunov Function Operations
 *===========================================================================*/

/**
 * @brief Evaluate a PWQ Lyapunov function at a point.
 *
 * V(x) = x^T P_i x + 2 q_i^T x + r_i for x in region_i
 *
 * @param lyap     PWQ Lyapunov function
 * @param x        State vector (n_state)
 * @param region   Region index containing x
 * @return V(x) value
 *
 * Complexity: O(n_state^2)
 */
double pwa_pwq_evaluate(const PWQLyapunov *lyap, const double *x, int region);

/**
 * @brief Evaluate the time derivative of a PWQ Lyapunov function.
 *
 * For CT-PWA in region i with dynamics dx/dt = A_i x + f_i:
 *   dV/dt = 2 x^T P_i (A_i x + f_i) + 2 q_i^T (A_i x + f_i)
 *
 * @param lyap     PWQ Lyapunov function
 * @param dyn      Affine dynamics for the current region
 * @param x        State vector
 * @param region   Region index
 * @return dV/dt value (negative ⇒ stability)
 */
double pwa_pwq_derivative(const PWQLyapunov *lyap,
                           const PWAAffineDynamics *dyn,
                           const double *x, int region);

/**
 * @brief Construct PWQ Lyapunov function from region matrices P_i, q_i, r_i.
 *
 * @param n_regions Number of regions
 * @param n_state   State dimension
 * @param P         Array of P_i matrices (each n_state × n_state)
 * @param q         Array of q_i vectors (each n_state)
 * @param r         Array of r_i scalars
 * @return New PWQ Lyapunov function, or NULL on error
 */
PWQLyapunov* pwa_pwq_create(int n_regions, int n_state,
                              const double **P, const double **q, const double *r);

/**
 * @brief Destroy a PWQ Lyapunov function.
 *
 * @param lyap PWQ Lyapunov function to destroy
 */
void pwa_pwq_destroy(PWQLyapunov *lyap);

/**
 * @brief Check positivity of PWQ Lyapunov function.
 *
 * Verifies V(x) > 0 for all x ∈ R_i \ {0}.
 * Uses S-procedure: ∃ U_i ≽ 0 s.t. P_i - H_i^T U_i H_i ≻ 0
 *
 * @param lyap PWQ Lyapunov function
 * @param sys  PWA system
 * @return 1 if positive definite, 0 otherwise
 */
int pwa_pwq_check_positive(const PWQLyapunov *lyap, const PWASystem *sys);

/**
 * @brief Check decrease condition for PWQ Lyapunov function.
 *
 * Verifies dV/dt < -αV along trajectories in each region.
 * Uses S-procedure approach.
 *
 * @param lyap  PWQ Lyapunov function
 * @param sys   PWA system
 * @param alpha Decay rate requirement
 * @return 1 if decreasing, 0 otherwise
 */
int pwa_pwq_check_decrease(const PWQLyapunov *lyap, const PWASystem *sys,
                            double alpha);

/**
 * @brief Check continuity of PWQ Lyapunov function across region boundaries.
 *
 * V is continuous at a boundary point x if for adjacent regions i,j:
 *   x^T(P_i - P_j)x + 2(q_i - q_j)^T x + (r_i - r_j) = 0
 *
 * @param lyap      PWQ Lyapunov function
 * @param sys       PWA system
 * @param tolerance Tolerance for continuity check
 * @return 1 if continuous, 0 otherwise
 */
int pwa_pwq_check_continuity(const PWQLyapunov *lyap, const PWASystem *sys,
                              double tolerance);

/*===========================================================================
 * L4: S-Procedure and LMI Methods
 *===========================================================================*/

/**
 * @brief Apply the S-procedure to a set of quadratic constraints.
 *
 * Given f_0(x) = x^T Q_0 x and f_k(x) = x^T Q_k x for k=1..m,
 * checks if f_0(x) ≥ 0 whenever f_k(x) ≥ 0 for all k.
 *
 * S-procedure sufficient condition:
 * ∃ τ_k ≥ 0 s.t. Q_0 - Σ τ_k Q_k ≽ 0 (positive semidefinite)
 *
 * @param Q0     Central quadratic form (n × n)
 * @param Qk     Constraint quadratic forms (m array of n × n)
 * @param n      Dimension
 * @param m      Number of constraints
 * @param tau    Output: S-procedure multipliers (length m)
 * @return 1 if S-procedure succeeds, 0 otherwise
 *
 * Reference: Yakubovich (1971), "S-procedure in nonlinear control theory."
 */
int pwa_s_procedure(const double *Q0, const double **Qk, int n, int m,
                     double *tau);

/**
 * @brief Build LMI for PWA stability analysis.
 *
 * Constructs the LMI: F_0 + Σ_vars x_i F_i ≽ 0
 * that encodes PWQ Lyapunov conditions for all regions.
 *
 * @param sys   PWA system
 * @param lmi   Output: LMI structure (caller allocates)
 * @return 0 on success, -1 on error
 */
int pwa_build_stability_lmi(const PWASystem *sys, PWALMI *lmi);

/**
 * @brief Check if a symmetric matrix is positive semidefinite.
 *
 * Uses Cholesky decomposition attempt.
 *
 * @param M Symmetric matrix (n × n, row-major)
 * @param n Dimension
 * @return 1 if M ≽ 0, 0 otherwise
 *
 * Complexity: O(n^3)
 */
int pwa_matrix_is_psd(const double *M, int n);

/**
 * @brief Check if a symmetric matrix is positive definite.
 *
 * @param M Symmetric matrix (n × n, row-major)
 * @param n Dimension
 * @return 1 if M ≻ 0, 0 otherwise
 */
int pwa_matrix_is_pd(const double *M, int n);

/*===========================================================================
 * L4: Invariant Set Computation
 *===========================================================================*/

/**
 * @brief Compute a polyhedral invariant set for a PWA system in one region.
 *
 * For DT: find polytope S s.t. A_i S + f_i ⊆ S
 * For CT: find polytope S s.t. for all x on boundary, vector field points inward
 *
 * @param sys   PWA system
 * @param d     Region index
 * @param iset  Output: invariant set (caller allocates poly)
 * @return 0 on success, -1 on no invariant set found
 *
 * Complexity: O(n_state^2 · n_vertices)
 */
int pwa_compute_invariant_set(const PWASystem *sys, int d,
                               PWAInvariantSet *iset);

/**
 * @brief Compute the maximal positively invariant set in a polyhedral region.
 *
 * Iterative algorithm: S_0 = region, S_{k+1} = Pre(S_k) ∩ S_k
 *
 * @param sys PWA system
 * @param d   Region index
 * @param max_iter Maximum iterations
 * @param iset Output: maximal invariant set
 * @return Number of iterations, -1 on error
 *
 * Complexity: O(max_iter · n_vertices^2)
 * Reference: Gilbert & Tan (1991), "Linear systems with state and
 *   control constraints: the maximal output admissible set."
 */
int pwa_maximal_invariant_set(const PWASystem *sys, int d, int max_iter,
                               PWAInvariantSet *iset);

/**
 * @brief Check if a polyhedron is a (controlled) invariant set.
 *
 * @param sys PWA system
 * @param d   Region index
 * @param S   Candidate invariant set
 * @return 1 if S is invariant, 0 otherwise
 */
int pwa_check_invariant(const PWASystem *sys, int d,
                         const PWAPolyhedron *S);

/*===========================================================================
 * L4: Common Lyapunov Function for Switched Systems
 *===========================================================================*/

/**
 * @brief Find a common quadratic Lyapunov function for all PWA regions.
 *
 * A common Lyapunov function V(x) = x^T P x (same P for all regions)
 * guarantees stability under arbitrary switching.
 *
 * Finds P ≻ 0 s.t. A_i^T P + P A_i ≺ 0 for all i (CT) or
 *                    A_i^T P A_i - P ≺ 0 for all i (DT).
 *
 * @param sys PWA system
 * @param P   Output: common P matrix (n_state × n_state, symmetric)
 * @return 1 if common Lyapunov function found, 0 otherwise
 *
 * Theorem (Liberzon 1999, p.34): A switched linear system is
 * GUAS under arbitrary switching iff all subsystems share a
 * common Lyapunov function.
 */
int pwa_common_lyapunov(const PWASystem *sys, double *P);

/**
 * @brief Check stability under minimum dwell time.
 *
 * Verifies that the switched system is stable if the time between
 * switches exceeds a minimum dwell time τ_d.
 *
 * @param sys       PWA system
 * @param tau_dwell Minimum dwell time
 * @return 1 if stable under dwell time, 0 otherwise
 *
 * Theorem (Morse 1996): If each subsystem is exponentially stable,
 * then ∃ τ_d > 0 s.t. the switched system is GUAS for dwell time ≥ τ_d.
 */
int pwa_dwell_time_stability(const PWASystem *sys, double tau_dwell);

#ifdef __cplusplus
}
#endif

#endif /* PWA_STABILITY_H */
