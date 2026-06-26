#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "switched_types.h"
#include "switched_applications.h"
#include "switched_stability.h"
#include "switched_dwell_time.h"

/**
 * Example 3: Networked Control with Packet Dropouts (L7)
 *
 * Demonstrates stability analysis of a networked control system
 * where packet dropouts cause mode switches between:
 *   - Mode 0: Packet received (feedback active, stable closed-loop)
 *   - Mode 1: Packet lost (open-loop, potentially unstable)
 *
 * Uses average dwell-time analysis to determine the maximum
 * allowable dropout rate for guaranteed stability.
 *
 * References:
 *   Zhang, W., Branicky, M.S. & Phillips, S.M. (2001).
 *   Stability of networked control systems. IEEE Control Systems.
 */
int main(void) {
    printf("=== Example 3: Networked Control with Packet Dropouts ===\n\n");
    srand((unsigned int)time(NULL));

    /* Plant: unstable first-order system, dx/dt = a*x + b*u
     * Controller: u = k*x (state feedback)
     * Closed-loop: dx/dt = (a + b*k)*x = a_cl * x (stable if a_cl < 0) */
    double a = 0.5;     /* Open-loop: unstable (eigenvalue = 0.5) */
    double b = 1.0;     /* Input gain */
    double k = -2.0;    /* Controller gain */
    double a_cl = a + b * k; /* Closed-loop eigenvalue = 0.5 + 1*(-2) = -1.5 */
    double dropout_rate = 0.2; /* 20% packet loss */

    printf("Plant: dx/dt = %.1f*x + %.1f*u\n", a, b);
    printf("Controller: u = %.1f*x\n", k);
    printf("Closed-loop: dx/dt = %.1f*x (stable)\n", a_cl);
    printf("Open-loop: dx/dt = %.1f*x (UNSTABLE)\n", a);
    printf("Packet dropout rate: %.0f%%\n\n", 100.0 * dropout_rate);

    /* Create switched system model: 2 modes */
    SwitchedSystem *sys = sss_create("NCS Dropout", 1, 2);

    /* Mode 0: Packet received, closed-loop stable */
    SwitchedMatrix A0 = sm_create(1, 1);
    sm_set(&A0, 0, 0, a_cl); /* -1.5, Hurwitz */
    sss_add_subsystem(sys, 0, &A0);

    /* Mode 1: Packet lost, open-loop unstable */
    SwitchedMatrix A1 = sm_create(1, 1);
    sm_set(&A1, 0, 0, a); /* 0.5, NOT Hurwitz */
    sss_add_subsystem(sys, 1, &A1);

    /* Switched stability analysis */
    printf("--- Switched System Analysis ---\n");

    /* Only Mode 0 is Hurwitz */
    bool gues = sss_is_gues_arbitrary(sys);
    printf("GUES under arbitrary switching: %s\n", gues ? "YES" : "NO");
    printf("(Mode 1 is unstable, so need dwell-time analysis)\n");

    /* Compute Lyapunov function for Mode 0 */
    /* A_cl^T P + P A_cl = -1 => P = 1/(2*|a_cl|) = 1/3 */
    double P_val = 1.0 / (2.0 * fabs(a_cl));
    printf("Lyapunov for Mode 0: P = %.4f, V(x) = %.4f*x^2\n", P_val, P_val);

    /* For Mode 1 (unstable), P is still positive definite:
     * A^T P + P A = 2*a*P = 2*0.5*0.333 = 0.333 > 0
     * V grows at rate 2*a = 1.0 */

    /* Dwell-time analysis:
     * lambda_0 = |a_cl| = 1.5 (decay in Mode 0)
     * mu = 1 (same P for both modes, but Mode 1 is unstable)
     * Need more sophisticated analysis */
    double lambda_0 = fabs(a_cl);
    double mu = 1.0;
    double tau_a_star = sdt_compute_avg_dwell(lambda_0, mu);
    printf("\nStability margin (Mode 0): lambda_0 = %.4f\n", lambda_0);
    printf("Required avg dwell time: tau_a* = %.4f s\n", tau_a_star);

    /* Probability-based stability:
     * If dropout rate = p, then fraction of time in Mode 0 = 1-p
     * Average dwell time in Mode 0: tau_0 = (1-p)/p * tau_1
     * For stability: tau_0 > tau_a* */
    double p_crit = 1.0 / (1.0 + tau_a_star);
    printf("Critical dropout rate for stability: p_crit = %.4f (%.1f%%)\n",
           p_crit, 100.0 * p_crit);

    /* NCS simulation */
    printf("\n--- NCS Simulation (10 seconds) ---\n");
    double ctrl_gain_arr[1] = {k};
    NetworkedControlDropout *ncs = ncs_create(1, ctrl_gain_arr, dropout_rate, (int)(1.0 / dropout_rate));
    ncs_simulate(ncs, 10.0, 0.01);

    /* MADB computation */
    SwitchedMatrix A_mat = sm_create(1, 1), B_mat = sm_create(1, 1), K_mat = sm_create(1, 1);
    sm_set(&A_mat, 0, 0, a);
    sm_set(&B_mat, 0, 0, b);
    sm_set(&K_mat, 0, 0, k);
    int madb = ncs_compute_madb(&A_mat, &B_mat, &K_mat, 1);
    printf("\nMADB (Max Allowable Dropout Bound): %d packets\n", madb);

    ncs_stability_analysis(ncs);

    /* Cleanup */
    sm_free(&A0); sm_free(&A1);
    sm_free(&A_mat); sm_free(&B_mat); sm_free(&K_mat);
    ncs_free(ncs);
    sss_free(sys);

    printf("\n=== Example 3 Complete ===\n");
    return 0;
}