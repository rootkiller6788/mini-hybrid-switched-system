/**
 * @file example_dcdc_converter_pwa.c
 * @brief DC-DC Buck Converter as PWA System — L6 Canonical Problem
 *
 * A DC-DC buck converter is a classic hybrid system: the switch
 * toggles between two modes (ON/OFF), making it naturally
 * piecewise affine. The continuous state is [i_L, v_C] (inductor
 * current, capacitor voltage).
 *
 * Mode ON:  L di/dt = -v_C + V_in,  C dv/dt = i_L - v_C/R
 * Mode OFF: L di/dt = -v_C,         C dv/dt = i_L - v_C/R
 *
 * Parameters: L=1mH, C=1mF, R=10Ω, V_in=10V
 *
 * Nine-school alignment:
 *   MIT 6.334 Power Electronics — Lec 3 Buck converter
 *   ETH 227-0220 Model Reduction — Lec 8 Hybrid models
 *   Berkeley EE222 — Ch 7 Switched systems
 */

#include "pwa_defs.h"
#include "pwa_simulation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(void)
{
    printf("=== DC-DC Buck Converter PWA Example ===\n\n");

    /* Parameters */
    double L = 1e-3;     /* 1 mH */
    double C = 1e-3;     /* 1 mF */
    double R = 10.0;     /* 10 Ohm */
    double Vin = 10.0;   /* 10 V */

    /* State: [i_L, v_C] */
    /* Mode ON:  di/dt = (-v_C + Vin)/L, dv/dt = (i_L - v_C/R)/C */
    /* Mode OFF: di/dt = -v_C/L,         dv/dt = (i_L - v_C/R)/C */

    PWASystem *sys = pwa_system_create(2, 0, 2, 2, 1, 0.0);
    if (!sys) { printf("Failed to create system\n"); return 1; }

    /* A_ON = [[0, -1/L], [1/C, -1/(RC)]] */
    double A_on[4] = {0.0, -1.0/L, 1.0/C, -1.0/(R*C)};
    double f_on[2] = {Vin / L, 0.0};

    /* A_OFF = [[0, -1/L], [1/C, -1/(RC)]] (same A, no Vin) */
    double A_off[4] = {0.0, -1.0/L, 1.0/C, -1.0/(R*C)};
    double f_off[2] = {0.0, 0.0};

    double C_out[4] = {1.0, 0.0, 0.0, 1.0};

    pwa_add_dynamics(sys, A_on, NULL, f_on, C_out, NULL, NULL);   /* Mode ON */
    pwa_add_dynamics(sys, A_off, NULL, f_off, C_out, NULL, NULL); /* Mode OFF */

    /* The switch is external (controlled). Regions defined by switch state.
     * For simulation, we use a PWM-like switching: region 0 when switch=ON,
     * region 1 when switch=OFF. The "input" is the duty cycle.
     *
     * For simplicity, define both regions covering the whole state space
     * and simulate with external mode switching. */

    /* Region 0 (ON): entire state space */
    double H0[8] = {1,0, -1,0, 0,1, 0,-1};
    double K0[4] = {100, 100, 100, 100};
    pwa_add_region(sys, H0, K0, 4, 0);

    /* Region 1 (OFF): entire state space */
    double H1[8] = {1,0, -1,0, 0,1, 0,-1};
    double K1[4] = {100, 100, 100, 100};
    pwa_add_region(sys, H1, K1, 4, 1);

    printf("Buck converter PWA model:\n");
    printf("  L=%.1e H, C=%.1e F, R=%.1f Ohm, Vin=%.1f V\n", L, C, R, Vin);
    pwa_system_print(sys, 1);

    /* Simulate switching sequence: PWM at 10 kHz, duty cycle = 0.5 */
    double T_sw = 1e-4;  /* 100 us switching period */
    double duty = 0.5;

    /* Use manual simulation with mode switching */
    double x0[2] = {0.0, 0.0};  /* Start from zero */
    double t = 0.0, t_end = 0.01;  /* 10 ms simulation */
    double dt = 1e-6;
    int max_steps = 20000;

    PWATrajectory *traj = pwa_trajectory_create(2, 0, 2, max_steps);
    if (!traj) { pwa_system_destroy(sys); return 1; }

    double x[2] = {x0[0], x0[1]};
    int cur_mode = 0;  /* Start ON */
    int step = 0;

    /* Store initial */
    traj->t_hist[0] = 0.0;
    traj->x_hist[0] = x[0]; traj->x_hist[1] = x[1];
    traj->region_hist[0] = cur_mode;
    traj->n_steps = 1;
    step = 1;

    while (t < t_end && step < max_steps) {
        /* Determine mode based on PWM */
        double t_in_period = fmod(t, T_sw);
        int desired_mode = (t_in_period < duty * T_sw) ? 0 : 1;

        if (desired_mode != cur_mode) {
            /* Mode switch event */
            cur_mode = desired_mode;
            if (traj->n_events < 100) {
                traj->events[traj->n_events].t = t;
                traj->events[traj->n_events].from_region = 1 - cur_mode;
                traj->events[traj->n_events].to_region = cur_mode;
                traj->n_events++;
            }
        }

        /* Integrate one step using RK4 */
        const PWAAffineDynamics *dyn = &sys->dynamics[cur_mode];
        /* Manual RK4 */
        double k1[2], k2[2], k3[2], k4[2], xt[2];

        /* k1 = f(x) */
        k1[0] = dyn->A[0]*x[0] + dyn->A[1]*x[1] + dyn->f[0];
        k1[1] = dyn->A[2]*x[0] + dyn->A[3]*x[1] + dyn->f[1];

        /* k2 */
        xt[0] = x[0] + 0.5*dt*k1[0]; xt[1] = x[1] + 0.5*dt*k1[1];
        k2[0] = dyn->A[0]*xt[0] + dyn->A[1]*xt[1] + dyn->f[0];
        k2[1] = dyn->A[2]*xt[0] + dyn->A[3]*xt[1] + dyn->f[1];

        /* k3 */
        xt[0] = x[0] + 0.5*dt*k2[0]; xt[1] = x[1] + 0.5*dt*k2[1];
        k3[0] = dyn->A[0]*xt[0] + dyn->A[1]*xt[1] + dyn->f[0];
        k3[1] = dyn->A[2]*xt[0] + dyn->A[3]*xt[1] + dyn->f[1];

        /* k4 */
        xt[0] = x[0] + dt*k3[0]; xt[1] = x[1] + dt*k3[1];
        k4[0] = dyn->A[0]*xt[0] + dyn->A[1]*xt[1] + dyn->f[0];
        k4[1] = dyn->A[2]*xt[0] + dyn->A[3]*xt[1] + dyn->f[1];

        x[0] += dt/6.0 * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        x[1] += dt/6.0 * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);

        t += dt;

        /* Store every 10th step */
        if (step % 10 == 0 && step < max_steps) {
            traj->t_hist[traj->n_steps] = t;
            traj->x_hist[traj->n_steps*2] = x[0];
            traj->x_hist[traj->n_steps*2+1] = x[1];
            traj->region_hist[traj->n_steps] = cur_mode;
            traj->n_steps++;
        }
        step++;
    }

    printf("Simulation: t=%.6f s, %d steps, v_C final=%.4f V, i_L final=%.4f A\n",
           t, traj->n_steps, x[1], x[0]);
    printf("Events: %d\n", traj->n_events);
    printf("Output voltage ripple: %.4f V (expected ~0.5V for 50%% duty)\n",
           x[1] - Vin * duty);

    pwa_trajectory_export_csv(traj, "build/dcdc_traj.csv");

    pwa_trajectory_destroy(traj);
    pwa_system_destroy(sys);
    printf("\nDone.\n");
    return 0;
}
