/**
 * @file example_bouncing_ball_pwa.c
 * @brief Bouncing Ball as PWA Hybrid System — L6 Canonical Problem
 *
 * The bouncing ball is the canonical example of a hybrid system
 * with state-dependent resets. When the ball hits the ground
 * (x1 = 0, x2 < 0), the velocity is reversed and scaled by
 * a restitution coefficient.
 *
 * State: x1 = height, x2 = velocity
 * Flow:  dx1/dt = x2, dx2/dt = -g
 * Guard: x1 ≤ 0 and x2 < 0
 * Reset: x2 → -c * x2 (where c is coefficient of restitution)
 *
 * This is a PWA system with one region (free fall) and a
 * state-dependent reset (jump) map at the boundary.
 *
 * Nine-school alignment:
 *   Berkeley EE222 — Ch 7 Hybrid systems
 *   Caltech CDS140 — Lec 16 Hybrid models
 *   MIT 6.832 Underactuated Robotics — Lec 5 Hybrid systems
 */

#include "pwa_defs.h"
#include "pwa_simulation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    printf("=== Bouncing Ball PWA Example ===\n\n");

    double g = 9.81;     /* Gravity */
    double c = 0.85;     /* Coefficient of restitution */
    double h0 = 10.0;    /* Initial height (m) */
    double t_end = 15.0; /* Simulation time */

    /* The system has one mode (free fall):
     * dx1/dt = x2
     * dx2/dt = -g
     *
     * With reset at boundary x1 ≤ 0: x2 → -c * x2
     *
     * We model this as a PWA system with 1 region covering the
     * state space x1 ≥ 0 (height non-negative). */

    PWASystem *sys = pwa_system_create(2, 0, 2, 1, 1, 0.0);
    if (!sys) { printf("Failed to create system\n"); return 1; }

    /* Free-fall dynamics */
    double A[4] = {0.0, 1.0, 0.0, 0.0};
    double f[2] = {0.0, -g};
    double C_out[4] = {1.0, 0.0, 0.0, 1.0};

    pwa_add_dynamics(sys, A, NULL, f, C_out, NULL, NULL);

    /* Region: x1 ≥ 0 → -x1 ≤ 0 */
    double H[4] = {-1, 0, 0, -1};  /* -x1 ≤ 0, -x2 ≤ 100 */
    double K[2] = {0, 100};
    pwa_add_region(sys, H, K, 2, 0);

    sys->x_min[0] = 0.0;
    sys->x_max[0] = 20.0;
    sys->x_min[1] = -30.0;
    sys->x_max[1] = 30.0;

    printf("Bouncing ball PWA model:\n");
    printf("  g=%.2f, restitution=%.2f, initial height=%.1f m\n", g, c, h0);
    pwa_system_print(sys, 1);

    /* Manual simulation with reset handling */
    PWASimConfig cfg = pwa_sim_config_default();
    cfg.t_end = t_end;
    cfg.dt_max = 0.005;
    cfg.max_steps = 20000;

    PWATrajectory *traj = pwa_trajectory_create(2, 0, 2, cfg.max_steps);
    if (!traj) { pwa_system_destroy(sys); return 1; }

    double t = 0.0;
    double x[2] = {h0, 0.0};
    double dt = cfg.dt_max;
    int step = 0;
    int bounce_count = 0;

    /* Store initial */
    traj->t_hist[0] = t;
    traj->x_hist[0] = x[0]; traj->x_hist[1] = x[1];
    traj->region_hist[0] = 0;
    traj->n_steps = 1;
    step = 1;

    while (t < t_end && step < cfg.max_steps) {
        double dt_use = dt;

        /* Detect if ground impact will occur within this step */
        /* Analytic: x1(t+dt) ≈ x1 + x2*dt - 0.5*g*dt^2 */
        double x1_pred = x[0] + x[1] * dt_use - 0.5 * g * dt_use * dt_use;

        if (x1_pred <= 0.0 && x[1] < 0.0) {
            /* Impact will occur. Find exact impact time.
             * Solve: x1 + v*t - 0.5*g*t^2 = 0 */
            double a = -0.5 * g;
            double b = x[1];
            double cc = x[0];
            double disc = b*b - 4*a*cc;
            if (disc >= 0) {
                double t_impact = (-b - sqrt(disc)) / (2*a);
                if (t_impact > 0 && t_impact < dt_use) {
                    dt_use = t_impact;
                }
            }
        }

        /* RK4 integration */
        double k1[2], k2[2], k3[2], k4[2];

        k1[0] = x[1];              k1[1] = -g;
        k2[0] = x[1] + 0.5*dt_use*k1[1];  k2[1] = -g;
        k3[0] = x[1] + 0.5*dt_use*k2[1];  k3[1] = -g;
        k4[0] = x[1] + dt_use*k3[1];      k4[1] = -g;

        x[0] += dt_use/6.0 * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        x[1] += dt_use/6.0 * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);

        t += dt_use;

        /* Check for impact */
        if (x[0] <= 1e-8 && x[1] < 0.0) {
            /* Reset: velocity reverses */
            x[1] = -c * x[1];
            x[0] = 0.0;
            bounce_count++;

            if (traj->n_events < 100) {
                traj->events[traj->n_events].t = t;
                traj->events[traj->n_events].from_region = 0;
                traj->events[traj->n_events].to_region = 0;
                traj->n_events++;
            }

            printf("  Bounce %d at t=%.4f s, v=%.4f m/s, peak=%.4f m\n",
                   bounce_count, t, fabs(x[1]),
                   x[1]*x[1]/(2*g));
        }

        if (step % 5 == 0) {
            traj->t_hist[traj->n_steps] = t;
            traj->x_hist[traj->n_steps*2] = x[0];
            traj->x_hist[traj->n_steps*2+1] = x[1];
            traj->region_hist[traj->n_steps] = 0;
            traj->n_steps++;
        }
        step++;
    }

    printf("\nSimulation complete: t=%.4f s, %d bounces\n", t, bounce_count);
    printf("Expected bounces (analytic): floor(log(h0/H_min)/log(1/c^2)) ≈ 6-8\n");

    /* Theoretical: total bounce time = sqrt(2h0/g) * (1+c)/(1-c) */
    double t_theory = sqrt(2*h0/g) * (1+c)/(1-c);
    printf("Theoretical total time: %.4f s\n", t_theory);

    pwa_trajectory_export_csv(traj, "build/bouncing_ball_traj.csv");
    printf("Trajectory exported to build/bouncing_ball_traj.csv\n");

    pwa_trajectory_destroy(traj);
    pwa_system_destroy(sys);
    printf("\nDone.\n");
    return 0;
}
