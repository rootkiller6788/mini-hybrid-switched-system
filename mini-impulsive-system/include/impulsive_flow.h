#ifndef IMPULSIVE_FLOW_H
#define IMPULSIVE_FLOW_H
/*
 * impulsive_flow.h — Continuous flow dynamics for impulsive systems
 *
 * Handles the continuous-time evolution between impulse instants:
 *   dx/dt = f(t, x),  for t ∈ (τ_k, τ_{k+1})
 *
 * References:
 *   Teschl (2012) "Ordinary Differential Equations and Dynamical Systems"
 *   Hairer, Nørsett, Wanner (1993) "Solving Ordinary Differential Equations I"
 */

#include "impulsive_types.h"

/* ── Linear Flow Systems ──────────────────────────────────────────────── */

/**
 * Linear flow: dx/dt = A * x
 *
 * FlowParams holds the A matrix (row-major) and its dimension.
 * The solution is x(t) = exp(A*(t-t0)) * x(t0).
 */
typedef struct {
    double *A;       /* [n*n] system matrix, row-major */
    int     n;       /* state dimension */
} ImpFlowLinear;

/**
 * affine_flow_create — Allocate linear flow parameters.
 *
 * @param A  [n*n] system matrix, copied internally
 * @param n  state dimension
 * @return allocated ImpFlowLinear or NULL on failure
 */
ImpFlowLinear* affine_flow_create(const double *A, int n);
void           affine_flow_free(ImpFlowLinear *flow);

/**
 * affine_flow_eval — Evaluate the vector field dx/dt = A*x.
 *
 * @param flow  linear flow descriptor
 * @param t     current time (ignored for autonomous systems)
 * @param x     state vector [flow->n]
 * @param dxdt  output: derivative [flow->n]
 * @param ctx   unused (pass NULL)
 * @return 0 on success
 */
int affine_flow_eval(double t, const double *x, int n, double *dxdt, void *ctx);

/* ── Nonlinear Flow Examples ──────────────────────────────────────────── */

/**
 * Nonlinear flow parameters for the Lorenz system:
 *   dx0/dt = sigma * (x1 - x0)
 *   dx1/dt = x0 * (rho - x2) - x1
 *   dx2/dt = x0 * x1 - beta * x2
 *
 * Classic chaotic system discovered by Edward Lorenz (1963).
 * Used as a testbed for impulsive synchronization.
 */
typedef struct {
    double sigma;  /* Prandtl number */
    double rho;    /* Rayleigh number */
    double beta;   /* geometric factor */
} ImpFlowLorenz;

ImpFlowLorenz* imp_flow_lorenz_create(double sigma, double rho, double beta);
void            imp_flow_lorenz_free(ImpFlowLorenz *lrz);

/**
 * imp_flow_lorenz_eval — Evaluate Lorenz vector field.
 *
 * Compatible with ImpVectorField signature; pass ImpFlowLorenz* as ctx.
 */
int imp_flow_lorenz_eval(double t, const double *x, int n,
                          double *dxdt, void *ctx);

/**
 * Nonlinear flow: Chua's circuit (canonical chaotic oscillator)
 *
 *   dx0/dt = alpha * (x1 - x0 - phi(x0))
 *   dx1/dt = x0 - x1 + x2
 *   dx2/dt = -beta * x1
 *   phi(x0) = m1*x0 + 0.5*(m0-m1)*(|x0+1| - |x0-1|)
 */
typedef struct {
    double alpha;
    double beta;
    double m0;
    double m1;
} ImpFlowChua;

ImpFlowChua* imp_flow_chua_create(double alpha, double beta, double m0, double m1);
void          imp_flow_chua_free(ImpFlowChua *chua);
int imp_flow_chua_eval(double t, const double *x, int n,
                        double *dxdt, void *ctx);

/**
 * Nonlinear flow: Duffing oscillator (forced)
 *   dx0/dt = x1
 *   dx1/dt = -delta*x1 + alpha*x0 - beta*x0^3 + gamma*cos(omega*t)
 */
typedef struct {
    double delta;
    double alpha;
    double beta;
    double gamma;
    double omega;
} ImpFlowDuffing;

ImpFlowDuffing* imp_flow_duffing_create(double delta, double alpha,
    double beta, double gamma, double omega);
void imp_flow_duffing_free(ImpFlowDuffing *duf);
int imp_flow_duffing_eval(double t, const double *x, int n,
                           double *dxdt, void *ctx);

/**
 * Nonlinear flow: Rössler attractor
 *   dx0/dt = -x1 - x2
 *   dx1/dt = x0 + a*x1
 *   dx2/dt = b + x2*(x0 - c)
 */
typedef struct {
    double a, b, c;
} ImpFlowRossler;

ImpFlowRossler* imp_flow_rossler_create(double a, double b, double c);
void imp_flow_rossler_free(ImpFlowRossler *r);
int imp_flow_rossler_eval(double t, const double *x, int n,
                           double *dxdt, void *ctx);

/**
 * Nonlinear flow: Van der Pol oscillator
 *   dx0/dt = x1
 *   dx1/dt = mu*(1 - x0^2)*x1 - x0
 */
typedef struct {
    double mu;
} ImpFlowVanDerPol;

ImpFlowVanDerPol* imp_flow_vdp_create(double mu);
void imp_flow_vdp_free(ImpFlowVanDerPol *vdp);
int imp_flow_vdp_eval(double t, const double *x, int n,
                       double *dxdt, void *ctx);

/**
 * Nonlinear flow: FitzHugh-Nagumo neuron model
 *   dv/dt = v - v^3/3 - w + I_ext
 *   dw/dt = epsilon*(v + a - b*w)
 */
typedef struct {
    double epsilon;
    double a;
    double b;
    double I_ext;
} ImpFlowFHN;

ImpFlowFHN* imp_flow_fhn_create(double epsilon, double a, double b, double I_ext);
void imp_flow_fhn_free(ImpFlowFHN *fhn);
int imp_flow_fhn_eval(double t, const double *x, int n,
                       double *dxdt, void *ctx);

/**
 * Nonlinear flow: Pendulum with damping
 *   dtheta/dt = omega
 *   domega/dt = -(g/L)*sin(theta) - b*omega
 */
typedef struct {
    double g_over_L;
    double b;
} ImpFlowPendulum;

ImpFlowPendulum* imp_flow_pendulum_create(double g_over_L, double b);
void imp_flow_pendulum_free(ImpFlowPendulum *pen);
int imp_flow_pendulum_eval(double t, const double *x, int n,
                            double *dxdt, void *ctx);

/**
 * Nonlinear flow: Simple harmonic oscillator
 *   dx0/dt = x1
 *   dx1/dt = -omega0^2 * x0
 */
typedef struct {
    double omega0_sq;
} ImpFlowHarmonic;

ImpFlowHarmonic* imp_flow_harmonic_create(double omega0);
void imp_flow_harmonic_free(ImpFlowHarmonic *h);
int imp_flow_harmonic_eval(double t, const double *x, int n,
                            double *dxdt, void *ctx);

/**
 * Generic flow wrapper for user-provided ImpVectorField with parameters.
 */
typedef struct {
    ImpVectorField  vf;
    void           *params;
    int             n;
} ImpFlowGeneric;

ImpFlowGeneric* imp_flow_generic_create(ImpVectorField vf, void *params, int n);
void imp_flow_generic_free(ImpFlowGeneric *gen);
int imp_flow_generic_eval(double t, const double *x, int n,
                           double *dxdt, void *ctx);

#endif /* IMPULSIVE_FLOW_H */
