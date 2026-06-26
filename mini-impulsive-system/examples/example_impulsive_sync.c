/* example_impulsive_sync.c -- Impulsive synchronization of two Lorenz systems.
 *
 * Master Lorenz system evolves autonomously; slave Lorenz system
 * receives impulsive corrections to synchronize.
 *
 * Master: dx/dt = Lorenz(x)
 * Slave:  dy/dt = Lorenz(y),  t != tau_k
 *         y^+ = y - B*(y - x), t = tau_k
 *
 * This demonstrates chaotic secure communication via impulsive sync.
 */
#include "impulsive_types.h"
#include "impulsive_flow.h"
#include "impulsive_jump.h"
#include "impulsive_sync.h"
#include "impulsive_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(void) {
    printf("=== Impulsive Synchronization Example ===\n\n");

    int n = 3;
    /* Lorenz parameters (chaotic regime) */
    double sigma = 10.0, rho = 28.0, beta = 8.0/3.0;
    ImpFlowLorenz *master_lorenz = imp_flow_lorenz_create(sigma, rho, beta);
    ImpFlowLorenz *slave_lorenz  = imp_flow_lorenz_create(sigma, rho, beta);

    /* Synchronization gain: B = 0.9*I (90% error correction per impulse) */
    double B[9] = {0.9, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0, 0.9};
    ImpSyncJumpLinear *sync_jump = imp_sync_jump_linear_create(B, n);

    /* Periodic impulses */
    double Delta = 0.1;
    ImpTimeSeq *seq = imp_time_seq_create_periodic(0.0, 20.0, Delta);

    /* Initial states set below in xm, xs arrays */

    ImpSolverConfig cfg = imp_solver_config_default();
    cfg.h_init = 0.001; cfg.method = IMP_SOLVER_RK4;

    double t = 0.0, h = cfg.h_init;
    double xm[3] = {1.0, 1.0, 1.0};
    double xs[3] = {5.0, -3.0, 10.0};
    int imp_idx = 0, step = 0;
    double max_err = 0.0, final_err = 0.0;

    while (t < seq->T && step < 200000) {
        double t_next = t + h;
        if (imp_idx < seq->count && t_next >= seq->times[imp_idx] - 1e-10) {
            h = seq->times[imp_idx] - t;
            /* Integrate master and slave */
            double dxm[3], dxs[3];
            imp_solver_rk4_step((ImpVectorField)imp_flow_lorenz_eval, master_lorenz,
                                t, xm, n, h, dxm);
            imp_solver_rk4_step((ImpVectorField)imp_flow_lorenz_eval, slave_lorenz,
                                t, xs, n, h, dxs);
            t += h;
            xm[0]=dxm[0]; xm[1]=dxm[1]; xm[2]=dxm[2];
            xs[0]=dxs[0]; xs[1]=dxs[1]; xs[2]=dxs[2];

            /* Apply sync impulse: correct slave error */
            double e[3], e_after[3];
            imp_sync_error_compute(xm, xs, n, IMP_SYNC_LINEAR_ERROR, e);
            imp_sync_jump_linear_eval(t, e, n, e_after, sync_jump);
            xs[0] = xm[0] + e_after[0];
            xs[1] = xm[1] + e_after[1];
            xs[2] = xm[2] + e_after[2];

            double err_norm = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
            if (err_norm > max_err) max_err = err_norm;
            final_err = err_norm;
            imp_idx++;
            h = cfg.h_init;

            if (step % 1000 == 0)
                printf("  t=%.3f  ||e||=%.6f\n", t, err_norm);
        } else {
            double dxm[3], dxs[3];
            imp_solver_rk4_step((ImpVectorField)imp_flow_lorenz_eval, master_lorenz,
                                t, xm, n, h, dxm);
            imp_solver_rk4_step((ImpVectorField)imp_flow_lorenz_eval, slave_lorenz,
                                t, xs, n, h, dxs);
            t += h;
            xm[0]=dxm[0]; xm[1]=dxm[1]; xm[2]=dxm[2];
            xs[0]=dxs[0]; xs[1]=dxs[1]; xs[2]=dxs[2];
        }
        step++;
    }

    printf("\nResults:\n");
    printf("  Steps: %d, Impulses: %d\n", step, imp_idx);
    printf("  Final slave:  [%.4f, %.4f, %.4f]\n", xs[0], xs[1], xs[2]);
    printf("  Final master: [%.4f, %.4f, %.4f]\n", xm[0], xm[1], xm[2]);
    printf("  Final ||e|| = %.6f\n", final_err);
    printf("  Max   ||e|| = %.6f\n", max_err);
    printf("  Synchronization: %s\n", (final_err < 0.1) ? "ACHIEVED" : "FAILED");

    imp_sync_jump_linear_free(sync_jump);
    imp_time_seq_free(seq);
    imp_flow_lorenz_free(master_lorenz);
    imp_flow_lorenz_free(slave_lorenz);
    return 0;
}
