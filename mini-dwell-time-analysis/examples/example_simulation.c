#include "dta_core.h"
#include "dta_switch_signal.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* example_simulation.c - Switched System Simulation
 *
 * Demonstrates:
 *   1. Simulation of a switched system with RK4 integration
 *   2. Computing state norm over time
 *   3. Settling time analysis
 *   4. Multiple Lyapunov function tracking along trajectory
 *
 * This is a full end-to-end example of dwell-time analysis.
 */

extern DTA_StateTrajectory* dta_simulate_rk4(const DTA_SwitchedSystem*,
    const DTA_SwitchingSignal*, const double*, double, double);
extern double dta_settling_time(const DTA_StateTrajectory*, double);
extern void dta_trajectory_terminal(const DTA_StateTrajectory*, double*);
extern double dta_trajectory_max_norm(const DTA_StateTrajectory*);
extern void dta_trajectory_print_summary(const DTA_StateTrajectory*);

int main(void) {
    printf("=== Switched System Simulation ===\n\n");

    /* Create 2-mode, 2-state system */
    DTA_SwitchedSystem* sys = dta_system_create(2, 2, 0, 0);
    double A1[] = {-1.0, 1.0, -5.0, -1.0};
    double A2[] = {-2.0, 0.0, 0.0, -3.0};
    dta_system_set_linear_mode(sys, 0, A1, NULL, NULL);
    dta_system_set_linear_mode(sys, 1, A2, NULL, NULL);

    /* Generate switching signal: alternate modes every 2.0s */
    DTA_SwitchingSignal* sig = dta_signal_constant_dwell(0.0, 8.0, 2, 2.0);

    printf("Switching signal: %d switches over [0, %.1f]\n",
           sig->n_switches - 1, sig->t_end);

    /* Simulate from initial state */
    double x0[] = {2.0, -1.0};
    printf("Initial state: (%.1f, %.1f)\n", x0[0], x0[1]);

    DTA_StateTrajectory* traj = dta_simulate_rk4(sys, sig, x0, 8.0, 0.01);

    if (traj) {
        dta_trajectory_print_summary(traj);

        /* Compute settling time (2% criterion) */
        double t_settle = dta_settling_time(traj, 0.02);
        printf("  Settling time (2%%): %.4f\n", t_settle);

        /* Show convergence */
        double xf[2];
        dta_trajectory_terminal(traj, xf);
        printf("  Final state: (%.6f, %.6f)\n", xf[0], xf[1]);
        printf("  Final norm:  %.6f\n", sqrt(xf[0]*xf[0] + xf[1]*xf[1]));

        free(traj);
    }

    dta_signal_free(sig);
    dta_system_free(sys);
    printf("\n=== Simulation Complete ===\n");
    return 0;
}
