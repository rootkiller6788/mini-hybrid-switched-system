#ifndef IMPULSIVE_SOLVER_H
#define IMPULSIVE_SOLVER_H
/*
 * impulsive_solver.h -- Numerical solvers for impulsive differential equations
 *
 * Solves:  dx/dt = f(t, x),  t != tau_k
 *          x(tau_k^+) = x(tau_k^-) + I_k(x(tau_k^-))
 *
 * Methods:
 *   - Forward Euler (fixed-step, 1st order)
 *   - RK4 (fixed-step, 4th order)
 *   - RK45 Dormand-Prince (adaptive step, 5th order)
 *   - Event detection via zero-crossing (for state-driven impulses)
 *
 * References:
 *   Hairer, Norsett, Wanner (1993) "Solving ODE I"
 *   Shampine & Gordon (1975) "Computer Solution of ODE"
 */

#include "impulsive_types.h"

/* -- Solver Configuration -- */

typedef enum {
    IMP_SOLVER_EULER    = 0,
    IMP_SOLVER_RK4      = 1,
    IMP_SOLVER_RK45     = 2,
    IMP_SOLVER_HEUN     = 3,
    IMP_SOLVER_MIDPOINT = 4
} ImpSolverMethod;

typedef struct {
    ImpSolverMethod method;
    double  h_init;
    double  h_min;
    double  h_max;
    double  atol;
    double  rtol;
    int     max_steps;
    bool    detect_events;
    double  event_tol;
} ImpSolverConfig;

ImpSolverConfig imp_solver_config_default(void);

/* -- Single-Step Integrators -- */

int imp_solver_euler_step(ImpVectorField f, void *ctx,
                           double t, const double *x, int n,
                           double h, double *x_next);

int imp_solver_rk4_step(ImpVectorField f, void *ctx,
                          double t, const double *x, int n,
                          double h, double *x_next);

int imp_solver_heun_step(ImpVectorField f, void *ctx,
                          double t, const double *x, int n,
                          double h, double *x_next);

int imp_solver_midpoint_step(ImpVectorField f, void *ctx,
                              double t, const double *x, int n,
                              double h, double *x_next);

/* -- Full Impulsive System Simulation -- */

int imp_solver_simulate(const ImpSystem *sys,
                         const double *x0, int n,
                         double t0, double T,
                         const ImpSolverConfig *cfg,
                         ImpSolution *sol);

/* -- Fixed-Step Simulation with Known Impulse Times -- */

int imp_solver_simulate_fixed(const ImpSystem *sys,
                               const double *x0, int n,
                               const ImpSolverConfig *cfg,
                               ImpSolution *sol);

/* -- Adaptive Step with Event Detection -- */

int imp_solver_simulate_adaptive(const ImpSystem *sys,
                                  const double *x0, int n,
                                  const ImpSolverConfig *cfg,
                                  ImpSolution *sol);

/* -- Dense Output / Interpolation -- */

int imp_solver_interpolate(const ImpSolution *sol, double t,
                            double *x_out, int n);

/* -- Error Estimation -- */

double imp_solver_estimate_error(const double *x1, const double *x2,
                                  int n, double atol, double rtol);

int imp_solver_adaptive_step(ImpVectorField f, void *ctx,
                              double t, const double *x, int n,
                              double *h, double *x_next,
                              double *err_out,
                              const ImpSolverConfig *cfg);

/* -- Apply Single Impulse -- */

int imp_solver_apply_impulse(ImpJumpMap I, void *ctx,
                              double t_k, double *x, int n);

#endif /* IMPULSIVE_SOLVER_H */
