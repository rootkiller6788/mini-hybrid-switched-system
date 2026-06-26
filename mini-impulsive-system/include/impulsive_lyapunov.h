#ifndef IMPULSIVE_LYAPUNOV_H
#define IMPULSIVE_LYAPUNOV_H
/*
 * impulsive_lyapunov.h -- Lyapunov stability theory for impulsive systems
 *
 * Lyapunov's direct method extended to impulsive systems.
 * A candidate Lyapunov function V: R^n -> R+ must satisfy:
 *
 *   1. V(0) = 0, V(x) > 0 for x != 0
 *   2. dV/dt <= -alpha * V  (during flow between impulses)
 *   3. V(x^+) <= rho * V(x^-)  (at impulse times)
 *
 * The dwell-time condition for stability:
 *   tau_D > ln(rho) / (-alpha)    if alpha < 0 (flow stable)
 *   tau_D < -ln(rho) / alpha       if alpha > 0 (need impulses)
 *
 * Key theorems:
 *   - Lyapunov stability for impulsive ODEs
 *   - Impulsive stabilization via Lyapunov control
 *   - Comparison principle for impulsive systems
 *   - Average dwell-time theorem (Hespanha & Morse, 1999)
 *
 * References:
 *   Yang (2001) "Impulsive Control Theory", Theorem 2.1-2.5
 *   Haddad et al. (2006) "Impulsive and Hybrid Dynamical Systems"
 *   Bainov & Simeonov (1989) "Systems with Impulse Effect"
 *   Lakshmikantham et al. (1989) "Theory of Impulsive Differential Equations"
 */

#include "impulsive_types.h"

struct ImpLyapunovFn {
    double (*V)(const double *x, int n, void *ctx);
    int    (*gradV)(const double *x, int n, double *grad, void *ctx);
    void    *ctx;
    char     name[64];
};

/* -- Quadratic Lyapunov Function -- */

typedef struct {
    double *P;
    int     n;
} ImpLyapQuad;

ImpLyapQuad* imp_lyap_quad_create(const double *P, int n);
void         imp_lyap_quad_free(ImpLyapQuad *lyap);
double       imp_lyap_quad_V(const double *x, int n, void *ctx);
int          imp_lyap_quad_gradV(const double *x, int n, double *grad, void *ctx);

double imp_lyap_quad_deriv_along_flow(const double *P, const double *A,
                                       const double *x, int n);
double imp_lyap_quad_jump_ratio(const double *P, const double *J,
                                 const double *x, int n);

/* -- Norm-Based Lyapunov Function -- */

typedef struct {
    int     p;
    double  q;
    int     n;
} ImpLyapNorm;

ImpLyapNorm* imp_lyap_norm_create(int p, double q, int n);
void         imp_lyap_norm_free(ImpLyapNorm *ln);
double       imp_lyap_norm_V(const double *x, int n, void *ctx);
int          imp_lyap_norm_gradV(const double *x, int n, double *grad, void *ctx);

/* -- Composite Lyapunov Function -- */

typedef struct {
    ImpLyapunovFn **components;
    double         *weights;
    int             k;
} ImpLyapComposite;

ImpLyapComposite* imp_lyap_composite_create(ImpLyapunovFn **comps,
                                              const double *w, int k);
void imp_lyap_composite_free(ImpLyapComposite *lc);
double imp_lyap_composite_V(const double *x, int n, void *ctx);
int imp_lyap_composite_gradV(const double *x, int n, double *grad, void *ctx);

/* -- Impulsive Lyapunov Stability Theorems -- */

int imp_lyap_check_stability(ImpLyapunovFn *lyap,
                              ImpVectorField f, void *f_ctx,
                              ImpJumpMap I, void *I_ctx,
                              const double *test_pts, int num_pts, int n,
                              double alpha, double rho,
                              ImpStabilityType *result);

double imp_lyap_dwell_time_bound(double c, double d);

double imp_lyap_estimate_contraction(ImpLyapunovFn *lyap,
                                      const double *x_before,
                                      const double *x_after, int n);

double imp_lyap_estimate_flow_rate(ImpLyapunovFn *lyap,
                                    const double *x_t,
                                    const double *x_t_dt,
                                    double dt, int n);

double imp_lyap_average_dwell_time(double lambda_val, double mu);

/* -- Lyapunov Matrix Equation Solvers -- */

int imp_lyap_solve_lyapunov_eq(const double *A, const double *Q,
                                double *P, int n);

int imp_lyap_solve_riccati_eq(const double *A, const double *B,
                               const double *Q, const double *R,
                               double *P, int n, int m);

int imp_lyap_check_pos_def(const double *P, int n);

#endif /* IMPULSIVE_LYAPUNOV_H */
