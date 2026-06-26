#ifndef DTA_MLF_H
#define DTA_MLF_H
#include "dta_core.h"

/* ==============================================================
 * dta_mlf.h - Multiple Lyapunov Function Construction
 *
 * Multiple Lyapunov functions (MLF) are the primary tool for
 * analyzing stability of switched systems.
 *
 * Key condition: For a family {V_i(x) = x^T P_i x}_{i=1..m},
 *   1. Each V_i decreases when mode i is active:
 *      V̇_i(x) = x^T (A_i^T P_i + P_i A_i) x < 0
 *   2. At switching instants t_k (σ changes from i to j):
 *      V_j(x(t_k)) ≤ V_i(x(t_k))
 *
 * The dwell time arises from the "jump" condition being relaxed:
 *   V_j(x(t_k)) ≤ μ V_i(x(t_k)),  μ ≥ 1
 * Then stability requires τ_d > ln(μ)/(2λ) where λ = min λ_i.
 *
 * References:
 *   Branicky (1998) "Multiple Lyapunov functions and other
 *     analysis tools for switched and hybrid systems", IEEE TAC
 *   DeCarlo, Branicky, Pettersson, Lennartson (2000)
 *     "Perspectives and results on the stability and
 *      stabilizability of hybrid systems", Proc. IEEE
 *   Liberzon (2003) Ch. 3.2 — MLF for stability analysis
 * ============================================================== */

/* --- MLF type --- */
typedef enum {
    DTA_MLF_QUADRATIC = 0,
    DTA_MLF_PIECEWISE_QUAD = 1,
    DTA_MLF_POLYHEDRAL = 2,
    DTA_MLF_SUM_OF_SQUARES = 3,
    DTA_MLF_PATH_DEPENDENT = 4
} DTA_MLFType;

/* --- Multiple Lyapunov function --- */
typedef struct {
    double** P_i;              /* P_i for each mode, each n×n */
    double* lambda_i;          /* Decay rate per mode */
    int n_modes;
    int n;                     /* State dimension */
    DTA_MLFType type;
    double mu;                 /* Coupling: max P_i eig / min P_j eig */
} DTA_MultipleLyapunov;

/* --- MLF verification result --- */
typedef struct {
    bool valid;
    bool all_pd;               /* All P_i positive definite */
    bool all_decrease;         /* V̇_i < 0 when mode active */
    bool switching_decrease;   /* V_j(t_k) ≤ μ V_i(t_k) at switches */
    double min_decay_rate;
    double max_coupling;
    double tau_d_required;     /* Minimum dwell time needed */
} DTA_MLFVerification;

/* --- API --- */

/** Allocate an MLF structure */
DTA_MultipleLyapunov* dta_mlf_create(int n_modes, int n);
void dta_mlf_free(DTA_MultipleLyapunov* mlf);

/** Set P_i for mode i */
int dta_mlf_set_P(DTA_MultipleLyapunov* mlf, int mode, const double* P);

/** Construct quadratic MLF by solving Lyapunov equations:
 *  A_i^T P_i + P_i A_i = -Q_i  for each mode.
 *  Q_i defaults to identity if NULL. */
DTA_MultipleLyapunov* dta_mlf_construct_quadratic(
    const DTA_SwitchedSystem* sys, const double** Q_matrices);

/** Construct MLF with optimized decay rates λ_i.
 *  Solves: P_i > 0, A_i^T P_i + P_i A_i + 2λ_i P_i < 0
 *  Maximizes min λ_i via iterative LMI. */
DTA_MultipleLyapunov* dta_mlf_construct_optimal(
    const DTA_SwitchedSystem* sys, int max_iter, double tol);

/** Compute the coupling constant μ = max_{i,j} λ_max(P_i P_j^{-1}) */
double dta_mlf_compute_mu(const DTA_MultipleLyapunov* mlf);

/** Evaluate V_i(x) = x^T P_i x */
double dta_mlf_evaluate(const DTA_MultipleLyapunov* mlf, int mode,
                         const double* x);

/** Evaluate V̇_i(x) = x^T (A_i^T P_i + P_i A_i) x */
double dta_mlf_derivative(const DTA_MultipleLyapunov* mlf, int mode,
                           const DTA_SwitchedSystem* sys, const double* x);

/** Verify that the MLF satisfies all conditions for GAS under dwell time τ_d */
DTA_MLFVerification dta_mlf_verify(const DTA_MultipleLyapunov* mlf,
    const DTA_SwitchedSystem* sys, const DTA_SwitchingSignal* sig,
    const double* x0, double tau_d);

/** Construct a common Lyapunov function V(x)=x^T P x if it exists.
 *  Checks: A_i^T P + P A_i < 0 for all i.
 *  Uses gradient descent on the cone of positive definite matrices. */
bool dta_mlf_construct_common(const DTA_SwitchedSystem* sys, double* P_out,
                               int max_iter, double step, double tol);

/** Piecewise quadratic MLF: V(x) = max_i { x^T P_i x } or min.
 *  Useful for state-dependent switching. */
DTA_MultipleLyapunov* dta_mlf_construct_piecewise(
    const DTA_SwitchedSystem* sys, double* region_centers,
    int n_regions);

/** Compute the dwell time required by a given MLF:
 *  τ_d* = ln(μ) / (2 min λ_i) */
double dta_mlf_required_dwell(const DTA_MultipleLyapunov* mlf);

/** Scale all P_i by factor: P_i ← α P_i. Changes mu by factor α. */
void dta_mlf_scale(DTA_MultipleLyapunov* mlf, double alpha);

#endif /* DTA_MLF_H */
