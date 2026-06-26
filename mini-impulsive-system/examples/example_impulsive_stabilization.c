/* example_impulsive_stabilization.c -- Impulsive stabilization of an unstable linear system.
 *
 * Demonstrates that an unstable continuous-time system can be stabilized
 * by applying impulsive control at discrete time instants.
 *
 * System: dx/dt = A*x (unstable, A has positive eigenvalue)
 * Impulses: x^+ = (I + J)*x  at periodic times tau_k = k * Delta
 *
 * The condition for stability: ||I+J|| < exp(-lambda_max(A) * Delta)
 */
#include "impulsive_types.h"
#include "impulsive_flow.h"
#include "impulsive_jump.h"
#include "impulsive_lyapunov.h"
#include "impulsive_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== Impulsive Stabilization Example ===\n\n");

    int n = 2;
    /* Unstable system: A has eigenvalues +0.5 (unstable) */
    double A[4] = {0.5, 0.0, 0.0, 0.5};
    ImpFlowLinear *flow = affine_flow_create(A, n);

    /* Stabilizing jump: J = -0.8*I reduces state by 80% */
    double J[4] = {-0.8, 0.0, 0.0, -0.8};
    ImpJumpLinear *jump = imp_jump_linear_create(J, n);

    /* Periodic impulses every Delta = 0.5 */
    double Delta = 0.5;
    ImpTimeSeq *seq = imp_time_seq_create_periodic(0.0, 5.0, Delta);
    ImpSystem *sys = imp_system_create(n,
        (ImpVectorField)affine_flow_eval,
        (ImpJumpMap)imp_jump_linear_eval, seq, NULL);
    sys->ctx = flow;
    /* Override I ctx via a wrapper -- for demo simplicity, use direct */
    sys->ctx = flow;
    sys->I = (ImpJumpMap)imp_jump_linear_eval;

    /* Simulate */
    double x0[2] = {10.0, 5.0};
    ImpSolverConfig cfg = imp_solver_config_default();
    cfg.method = IMP_SOLVER_RK4;
    cfg.h_init = 0.001;

    ImpSolution *sol = imp_solution_create(10000, 100, n);

    /* Manual simulation to handle jump context */
    double t = seq->t0, h = cfg.h_init;
    double *x = (double*)malloc((size_t)n * sizeof(double));
    x[0] = x0[0]; x[1] = x0[1];
    imp_solution_add_point(sol, t, x);

    int imp_idx = 0, step = 0;
    while (t < seq->T && step < 50000) {
        double t_next = t + h;
        if (imp_idx < seq->count && t_next >= seq->times[imp_idx] - 1e-10) {
            h = seq->times[imp_idx] - t;
            double xn[2];
            imp_solver_rk4_step(sys->f, flow, t, x, n, h, xn);
            t += h;
            x[0] = xn[0]; x[1] = xn[1];

            /* Apply jump */
            double x_before[2] = {x[0], x[1]};
            imp_jump_linear_eval(t, x_before, n, x, jump);
            imp_solution_add_jump(sol, x_before, x);
            imp_idx++;
            h = cfg.h_init;
        } else {
            double xn[2];
            imp_solver_rk4_step(sys->f, flow, t, x, n, h, xn);
            t += h;
            x[0] = xn[0]; x[1] = xn[1];
        }
        imp_solution_add_point(sol, t, x);
        step++;
    }

    printf("Simulation: %d steps, %d impulses\n", step, imp_idx);
    printf("  t=0:     x = [%.4f, %.4f]\n", x0[0], x0[1]);
    printf("  t=%.1f:  x = [%.4f, %.4f]\n", sol->T,
           sol->x[(size_t)(sol->npts-1) * n], sol->x[(size_t)(sol->npts-1) * n + 1]);

    double norm_final = sqrt(x[0]*x[0] + x[1]*x[1]);
    printf("  ||x(T)|| = %.6f\n", norm_final);
    printf("  Stabilization: %s\n", (norm_final < 1.0) ? "SUCCESS" : "FAILED");

    free(x);
    imp_solution_free(sol);
    imp_system_free(sys);
    imp_time_seq_free(seq);
    affine_flow_free(flow);
    imp_jump_linear_free(jump);
    return 0;
}
