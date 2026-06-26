/**
 * @file example_saturation_pwa.c
 * @brief PWA Model of Actuator Saturation — L6 Canonical Problem
 *
 * Actuator saturation is the simplest and most common nonlinearity
 * in control systems. The saturation function sat(u) = clip(u, -u_max, u_max)
 * is naturally piecewise affine with three regions.
 *
 * This example constructs a PWA model of a double integrator with
 * saturated input, simulates it, and demonstrates that the PWA
 * framework can exactly represent this hybrid behavior.
 *
 * System:  x1_dot = x2
 *          x2_dot = sat(u),   where sat(u) = { -1 if u < -1
 *                                              u  if -1 ≤ u ≤ 1
 *                                               1 if u > 1 }
 *
 * Nine-school alignment:
 *   MIT 6.241J — Lec 5 Saturation
 *   Stanford AA203 — Lec 8 Input constraints
 *   Berkeley EE222 — Ch 3 Saturation
 */

#include "pwa_defs.h"
#include "pwa_simulation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void)
{
    printf("=== PWA Saturation Example ===\n\n");

    /* Create a 2-state, 1-input, 2-output PWA system with 3 regions.
     * Regions correspond to:
     *   Region 0: u < -1  →  sat(u) = -1
     *   Region 1: -1 ≤ u ≤ 1  →  sat(u) = u
     *   Region 2: u > 1  →  sat(u) = 1
     *
     * State x = [position, velocity], Input u = commanded force
     * dx1/dt = x2
     * dx2/dt = sat(u)
     *
     * In matrix form: A = [[0,1],[0,0]], B = [[0],[0]], f = [0, sat(u)]
     * We encode sat(u) via the affine offset f_i:
     *   Region 0: f = [0, -1]
     *   Region 1: f = [0, 0], B = [[0],[1]]  →  dx2/dt = u
     *   Region 2: f = [0, 1]
     */

    PWASystem *sys = pwa_system_create(2, 1, 2, 3, 1, 0.0);
    if (!sys) { printf("Failed to create system\n"); return 1; }

    /* A = [[0,1],[0,0]] for all regions */
    double A[4] = {0.0, 1.0, 0.0, 0.0};
    double B0[2] = {0.0, 0.0};
    double B1[2] = {0.0, 1.0};
    double f0[2] = {0.0, -1.0};   /* Region 0: sat = -1 */
    double f1[2] = {0.0, 0.0};    /* Region 1: sat = u */
    double f2[2] = {0.0, 1.0};    /* Region 2: sat = 1 */
    double C[4] = {1.0, 0.0, 0.0, 1.0};

    /* Add dynamics for each region */
    pwa_add_dynamics(sys, A, B0, f0, C, NULL, NULL); /* Region 0 */
    pwa_add_dynamics(sys, A, B1, f1, C, NULL, NULL); /* Region 1 */
    pwa_add_dynamics(sys, A, B0, f2, C, NULL, NULL); /* Region 2 */

    /* Region constraints on [x1, x2, u]:
     * Region 0: u ≤ -1
     * Region 1: -1 ≤ u ≤ 1  (i.e., u ≤ 1, -u ≤ 1)
     * Region 2: u ≥ 1  (i.e., -u ≤ -1) */

    int nz = 3;
    /* Region 0: u ≤ -1 →  [0,0,1]·z ≤ -1 */
    double H0[3] = {0, 0, 1}; double K0[1] = {-1.0};
    pwa_add_region(sys, H0, K0, 1, 0);

    /* Region 1: -1 ≤ u ≤ 1 → u ≤ 1, -u ≤ 1 */
    double H1[6] = {0,0,1, 0,0,-1}; double K1[2] = {1.0, 1.0};
    pwa_add_region(sys, H1, K1, 2, 1);

    /* Region 2: u ≥ 1 → -u ≤ -1 */
    double H2[3] = {0,0,-1}; double K2[1] = {-1.0};
    pwa_add_region(sys, H2, K2, 1, 2);

    /* Set domain bounds */
    for (int i = 0; i < 2; i++) {
        sys->x_min[i] = -5.0;
        sys->x_max[i] = 5.0;
    }
    sys->u_min[0] = -3.0;
    sys->u_max[0] = 3.0;

    printf("PWA system with saturation created:\n");
    pwa_system_print(sys, 1);

    /* Simulate with sinusoidal input */
    printf("\nSimulating with u(t) = 2*sin(t) ...\n");
    double x0[2] = {1.0, 0.0};

    PWASimConfig cfg = pwa_sim_config_default();
    cfg.t_end = 20.0;
    cfg.dt_max = 0.02;
    cfg.max_steps = 2000;

    PWATrajectory *traj = pwa_trajectory_create(2, 1, 2, cfg.max_steps);
    if (!traj) { pwa_system_destroy(sys); return 1; }

    pwa_simulate_ct(sys, x0, NULL, NULL, &cfg, traj);

    printf("Simulation complete. %d steps, %d events\n",
           traj->n_steps, traj->n_events);
    pwa_trajectory_print(traj);

    /* Export to CSV */
    pwa_trajectory_export_csv(traj, "build/saturation_traj.csv");
    printf("Trajectory exported to build/saturation_traj.csv\n");

    /* Check for saturation regions visited */
    int visited[3] = {0, 0, 0};
    for (int i = 0; i < traj->n_steps; i++) {
        int r = traj->region_hist[i];
        if (r >= 0 && r < 3) visited[r] = 1;
    }
    printf("Regions visited: neg_sat=%d linear=%d pos_sat=%d\n",
           visited[0], visited[1], visited[2]);

    pwa_trajectory_destroy(traj);
    pwa_system_destroy(sys);
    printf("\nDone.\n");
    return 0;
}
