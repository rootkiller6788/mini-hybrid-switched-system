#ifndef IMPULSIVE_JUMP_H
#define IMPULSIVE_JUMP_H
/*
 * impulsive_jump.h — Impulsive reset / jump map definitions and operations
 *
 * Jump maps describe the instantaneous state change at impulse instants τ_k:
 *   x(τ_k^+) = x(τ_k^-) + I_k(x(τ_k^-))
 *
 * or equivalently:  Δx|_τ_k = I_k(x(τ_k^-))
 *
 * References:
 *   Bainov & Simeonov (1989) "Systems with Impulse Effect"
 *   Lakshmikantham et al. (1989) "Theory of Impulsive Differential Equations"
 *   Yang (2001) "Impulsive Control Theory", Springer
 */

#include "impulsive_types.h"

/* ── Linear Jump Map Parameters ───────────────────────────────────────── */

/**
 * ImpJumpLinear — Linear reset: Δx = J * x  (x^+ = (I+J) * x^-)
 *
 * J is stored as [n*n] row-major matrix.
 * After jump: x^+ = x^- + J*x^- = (I + J)*x^-
 */
typedef struct {
    double *J;       /* [n*n] jump matrix */
    int     n;       /* state dimension */
} ImpJumpLinear;

ImpJumpLinear* imp_jump_linear_create(const double *J, int n);
void           imp_jump_linear_free(ImpJumpLinear *jl);
int imp_jump_linear_eval(double t_k, const double *x_before, int n,
                          double *x_after, void *ctx);

/* ── Affine Jump Map Parameters ───────────────────────────────────────── */

/**
 * ImpJumpAffine — Affine reset: Δx = J*x + d  (x^+ = (I+J)*x + d)
 */
typedef struct {
    double *J;       /* [n*n] matrix */
    double *d;       /* [n] offset vector */
    int     n;
} ImpJumpAffine;

ImpJumpAffine* imp_jump_affine_create(const double *J, const double *d, int n);
void           imp_jump_affine_free(ImpJumpAffine *ja);
int imp_jump_affine_eval(double t_k, const double *x_before, int n,
                          double *x_after, void *ctx);

/* ── Projection Jump Map ──────────────────────────────────────────────── */

/**
 * ImpJumpProject — Euclidean projection onto a convex set C
 *
 * x^+ = argmin_{z in C} ||z - x^-||_2
 *
 * Currently supports: box constraints [lb[i], ub[i]]
 */
typedef struct {
    double *lb;      /* [n] lower bounds */
    double *ub;      /* [n] upper bounds */
    int     n;
} ImpJumpProject;

ImpJumpProject* imp_jump_project_create(const double *lb, const double *ub, int n);
void            imp_jump_project_free(ImpJumpProject *jp);
int imp_jump_project_eval(double t_k, const double *x_before, int n,
                           double *x_after, void *ctx);

/* ── Hard Reset Jump Map ──────────────────────────────────────────────── */

/**
 * ImpJumpHardReset — Fixed-target reset: x^+ = x_target (regardless of x^-)
 */
typedef struct {
    double *x_target;  /* [n] target state */
    int     n;
} ImpJumpHardReset;

ImpJumpHardReset* imp_jump_hard_reset_create(const double *x_target, int n);
void              imp_jump_hard_reset_free(ImpJumpHardReset *hr);
int imp_jump_hard_reset_eval(double t_k, const double *x_before, int n,
                              double *x_after, void *ctx);

/* ── Impulsive Control Map ────────────────────────────────────────────── */

/**
 * ImpJumpControl — Impulsive control: Δx = B_k * u_k
 *
 * The control input u_k is applied as an impulse at time τ_k,
 * causing an instantaneous state change.
 *
 * x^+ = x^- + B*u
 */
typedef struct {
    double *B;       /* [n*m] input matrix */
    double *u;       /* [m] control input vector (can be updated) */
    int     n;       /* state dimension */
    int     m;       /* input dimension */
} ImpJumpControl;

ImpJumpControl* imp_jump_control_create(const double *B, int n, int m);
void            imp_jump_control_free(ImpJumpControl *jc);
void            imp_jump_control_set_u(ImpJumpControl *jc, const double *u_new);
int imp_jump_control_eval(double t_k, const double *x_before, int n,
                           double *x_after, void *ctx);

/* ── State-Dependent Nonlinear Jump ───────────────────────────────────── */

/**
 * ImpJumpNonlinear — User-supplied nonlinear jump: Δx = g(x)
 *
 * The user provides a callback of type ImpJumpMap.
 */
typedef struct {
    ImpJumpMap  g;       /* nonlinear map x_before -> x_after */
    void       *params;  /* user parameters */
    int         n;
} ImpJumpNonlinear;

ImpJumpNonlinear* imp_jump_nonlinear_create(ImpJumpMap g, void *params, int n);
void              imp_jump_nonlinear_free(ImpJumpNonlinear *jn);
int imp_jump_nonlinear_eval(double t_k, const double *x_before, int n,
                             double *x_after, void *ctx);

/* ── Time-Varying Jump ────────────────────────────────────────────────── */

/**
 * ImpJumpTimeVarying — Jump map that depends on the impulse time t_k
 *
 * Useful for periodic or scheduled impulse control.
 */
typedef struct {
    double *J_k;     /* [n*n] time-dependent gain (updated per impulse) */
    int     n;
} ImpJumpTimeVarying;

ImpJumpTimeVarying* imp_jump_time_varying_create(int n);
void                imp_jump_time_varying_free(ImpJumpTimeVarying *jtv);
void                imp_jump_time_varying_set_gain(ImpJumpTimeVarying *jtv,
                                                    const double *J_k);
int imp_jump_time_varying_eval(double t_k, const double *x_before, int n,
                                double *x_after, void *ctx);

/* ── Jump Map Composition ─────────────────────────────────────────────── */

/**
 * imp_jump_compose — Compose two jump maps: x^+ = J2(J1(x^-))
 *
 * Applies J1 then J2 at the same impulse instant.
 * @return 0 on success
 */
int imp_jump_compose(ImpJumpMap J1, void *ctx1,
                      ImpJumpMap J2, void *ctx2,
                      double t_k, const double *x_before, int n,
                      double *x_after);

/**
 * imp_jump_sequence — Apply a sequence of jump maps at successive times.
 *
 * @param J       array of jump maps
 * @param ctx_arr array of contexts
 * @param num     number of jump maps
 */
int imp_jump_sequence(ImpJumpMap *J, void **ctx_arr, int num,
                       double t_k, const double *x_before, int n,
                       double *x_after);

#endif /* IMPULSIVE_JUMP_H */
