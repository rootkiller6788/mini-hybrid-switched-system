#ifndef DTA_LMI_H
#define DTA_LMI_H
#include "dta_core.h"

/* ==============================================================
 * dta_lmi.h - LMI-Based Dwell Time Computation
 *
 * Linear Matrix Inequalities (LMIs) are the computational engine
 * for dwell-time analysis. The key LMI for dwell time τ_d is:
 *
 * Find P_i = P_i^T > 0,  i = 1..m  such that:
 *   A_i^T P_i + P_i A_i + 2λ P_i < 0     (mode decay)
 *   P_i ≤ μ P_j,  μ = e^{2λ τ_d}          (coupling)
 *
 * These are solved via semidefinite programming (SDP).
 * This module provides a self-contained interior-point
 * implementation for small-scale problems.
 *
 * References:
 *   Boyd, El Ghaoui, Feron, Balakrishnan (1994)
 *     "Linear Matrix Inequalities in System and Control Theory"
 *   Geromel & Colaneri (2006) SIAM J. Control Optim.
 *   Nesterov & Nemirovskii (1994) "Interior-Point Polynomial
 *     Algorithms in Convex Programming"
 * ============================================================== */

/* --- LMI constraint type --- */
typedef enum {
    DTA_LMI_LYAPUNOV = 0,      /* A^T P + P A + Q < 0 */
    DTA_LMI_COUPLING = 1,      /* P_i - μ P_j ≤ 0 */
    DTA_LMI_POS_DEF = 2,       /* P > 0 */
    DTA_LMI_DWELL_TIME = 3     /* Full dwell-time LMI system */
} DTA_LMI_Type;

/* --- LMI problem description --- */
typedef struct {
    int n_vars;                /* Number of decision variables */
    int n_constraints;         /* Number of LMI constraints */
    int* constraint_dims;      /* Dimension of each LMI, length n_constraints */
    double** F_matrices;       /* Coefficient matrices, row-major */
    double* c;                 /* Objective: min c^T x */
    DTA_LMI_Type* constraint_types;
} DTA_LMI_Problem;

/* --- LMI solution --- */
typedef struct {
    double* x;                 /* Primal solution */
    double* Z;                 /* Dual solution (slack) matrices */
    double obj_value;          /* Optimal objective value */
    bool feasible;             /* Is the LMI feasible? */
    int iterations;            /* Interior-point iterations */
    double duality_gap;
    int n_vars;
} DTA_LMI_Solution;

/* --- Dwell-time specific LMI structures --- */
typedef struct {
    double** P_i;              /* Lyapunov matrices P_i, each n×n */
    double* lambda_i;          /* Decay rates per mode */
    double mu;                 /* Coupling constant */
    bool feasible;
    double tau_d;              /* Dwell time achieved */
} DTA_DwellLMI_Result;

/* --- API --- */

/** Create an LMI problem */
DTA_LMI_Problem* dta_lmi_problem_create(int n_vars, int n_constraints);
void dta_lmi_problem_free(DTA_LMI_Problem* prob);

/** Add a Lyapunov LMI: A^T P + P A + 2λ P < 0 */
int dta_lmi_add_lyapunov_constraint(DTA_LMI_Problem* prob,
    const double* A, int n, double lambda, int P_start_idx);

/** Add a coupling constraint: P_i - μ P_j ≤ 0 */
int dta_lmi_add_coupling_constraint(DTA_LMI_Problem* prob,
    int P_i_start, int P_j_start, int n, double mu);

/** Add P > 0 constraint (positive definiteness) */
int dta_lmi_add_posdef_constraint(DTA_LMI_Problem* prob,
    int P_start_idx, int n, double eps);

/** Solve LMI using primal-dual interior-point method.
 *  Returns solution with x and feasibility flag. */
DTA_LMI_Solution dta_lmi_solve(const DTA_LMI_Problem* prob,
                                int max_iter, double tol);

/** Check dwell time feasibility: does there exist MLF for given τ_d?
 *  This is the core dwell-time feasibility check. */
bool dta_lmi_check_dwell_feasibility(const DTA_SwitchedSystem* sys,
                                      double tau_d);

/** Solve the full dwell-time LMI for a switched linear system.
 *  Returns P_i matrices, λ_i, μ in the result. */
DTA_DwellLMI_Result dta_lmi_solve_dwell_time(
    const DTA_SwitchedSystem* sys, double tau_d);

/** Minimum dwell time via bisection + LMI feasibility:
 *  τ_d* = inf { τ_d > 0 : LMI feasible } */
double dta_lmi_min_dwell_bisection(const DTA_SwitchedSystem* sys,
    double tau_low, double tau_high, double tol, int max_iter);

/** Solve for common Lyapunov function via LMI:
 *  Find P > 0 s.t. A_i^T P + P A_i < 0 for all i */
bool dta_lmi_common_lyapunov(const DTA_SwitchedSystem* sys,
                              double* P_out);

/** Solve the Lyapunov-Metzler inequalities:
 *  Find P_i > 0 and Metzler matrix Π = (π_ij) with π_ij ≥ 0,
 *  Σ_j π_ij = 0, such that:
 *  A_i^T P_i + P_i A_i + Σ_j π_ji P_j < 0  for all i */
bool dta_lmi_lyapunov_metzler(const DTA_SwitchedSystem* sys,
    double** P_out, double** Pi_out);

/** Check positive definiteness of a symmetric matrix.
 *  Uses Cholesky decomposition attempt. */
bool dta_lmi_is_positive_definite(const double* M, int n, double tol);

/** Check negative definiteness: all eigenvalues < -tol */
bool dta_lmi_is_negative_definite(const double* M, int n, double tol);

/** Eigenvalue-based feasibility check for small systems (n≤4).
 *  Direct algebraic check without full SDP solver. */
bool dta_lmi_feasibility_small(const DTA_SwitchedSystem* sys,
                                double tau_d);

#endif /* DTA_LMI_H */
