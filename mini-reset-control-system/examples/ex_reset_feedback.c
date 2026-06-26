/* ex_reset_feedback.c - Reset Feedback Control Demo
 *
 * Demonstrates a complete reset feedback control system: plant +
 * reset controller in unity feedback configuration. Shows the
 * hybrid simulation of continuous flow + discrete reset jumps.
 *
 * Plant: first-order stable system G(s) = 1/(s + 1)
 * Controller: Clegg integrator
 *
 * This example shows:
 *   - Closed-loop step response with reset control
 *   - Hybrid simulation (flow + jump)
 *   - Reset event timing and Zeno detection
 *
 * Ref: Banos & Barreiro (2012) Chapter 4; Nesic et al. (2008)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "reset_core.h"
#include "reset_system.h"
#include "reset_simulation.h"

int main(void)
{
    printf("=== Reset Feedback Control Demo ===\n");


    /* Build plant: G(s) = 1/(s+1), state-space: A=-1, B=1, C=1, D=0 */
    ResetLinearBase *plant = reset_base_create(1, 1, 1);
    if (!plant) {
        printf("ERROR: Failed to create plant\n");

        return 1;
    }
    plant->A[0] = -1.0;
    plant->B[0] = 1.0;
    plant->C[0] = 1.0;
    /* D[0] = 0.0 (already zero) */

    printf("Plant: G(s) = 1/(s+1) (stable first-order)\n");

    printf("Plant Hurwitz stable: %s\n",

reset_is_hurwitz(plant) ? "YES" : "NO");


    /* Build reset controller: Clegg integrator */
    ResetSystem *ctrl = reset_clegg_create();
    if (!ctrl) {
        printf("ERROR: Failed to create controller\n");

        reset_base_free(plant);
        return 1;
    }
    printf("Controller: Clegg integrator (full reset on ZC)\n");

    printf("Controller reset ratio: %.2f\n", ctrl->ratio->rho);


    /* Create feedback loop */
    ResetFeedbackLoop *loop = reset_feedback_create(ctrl, plant);
    if (!loop) {
        printf("ERROR: Failed to create feedback loop\n");

        reset_sys_free(ctrl);
        reset_base_free(plant);
        return 1;
    }

    /* Simulation parameters */
    double dt = 0.001;
    double t_end = 10.0;
    int n_steps = (int)(t_end / dt);
    double *y_out = (double*)malloc((size_t)n_steps * sizeof(double));
    double *t_span = (double*)malloc((size_t)n_steps * sizeof(double));
    double *r_vals = (double*)malloc((size_t)n_steps * sizeof(double));

    if (!y_out || !t_span || !r_vals) {
        printf("ERROR: Memory allocation failed\n");

        return 1;
    }

    /* Unit step reference */
    for (int k = 0; k < n_steps; k++) {
        t_span[k] = k * dt;
        r_vals[k] = 1.0;
    }

    printf("Simulating step response (dt=%.4f, T=%.1f s, %d steps)...\n",

dt, t_end, n_steps);

    printf("Time(s)	Output	Error\n");

    printf("-------	------	-----\n");


    int n_resets = reset_feedback_simulate(loop, t_span, r_vals, y_out, n_steps);

    /* Print selected data points */
    int print_interval = (int)(0.5 / dt);
    for (int k = 0; k < n_steps; k += print_interval) {
        printf("%.3f	%.4f	%.4f\n",

t_span[k], y_out[k], 1.0 - y_out[k]);

    }

    /* Stability check */
    bool base_stable = reset_check_hbeta_stability(plant, ctrl);
    printf("Base closed-loop Hurwitz stable: %s\n",

base_stable ? "YES" : "NO");


    /* H-infinity norm */
    double hinf = reset_hinf_norm(plant, 50, 1e-6);
    printf("Plant H-infinity norm: %.4f\n", hinf);


    /* Frequency response at key frequency */
    double w_cross = 1.0;
    double sens = reset_sensitivity(plant, ctrl, w_cross);
    double comp = reset_complementary_sensitivity(plant, ctrl, w_cross);
    printf("At w=%.1f: Sensitivity=%.4f, Comp.Sensitivity=%.4f\n",

w_cross, sens, comp);


    printf("Total simulation steps: %d\n", n_steps);

    printf("Total reset events during simulation: %d\n", n_resets);


    const ResetInterval *ri = reset_get_interval_stats(ctrl);
    if (ri && ri->n_resets > 0) {
        printf("Reset interval stats:\n");

        printf("  Min interval: %.6f s\n", ri->t_min_interval);

        printf("  Max interval: %.6f s\n", ri->t_max_interval);

        printf("  Avg interval: %.6f s\n", ri->t_avg_interval);

    }

    printf("=== Demo Complete ===\n");


    free(y_out); free(t_span); free(r_vals);
    reset_feedback_free(loop);
    reset_sys_free(ctrl);
    reset_base_free(plant);

    return 0;
}
