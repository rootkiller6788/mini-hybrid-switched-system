#include "etc_core.h"
#include "etc_trigger.h"
#include "etc_dynamics.h"
#include <math.h>
#include <stdio.h>

/* ============================================================================
 * Example 1: First-Order Stabilization with Event-Triggered Control
 *
 * Plant:  ẋ = a x + b u   (first-order scalar system)
 * Controller: u = −k x   (proportional state feedback)
 *
 * With a = 1 (unstable), b = 1, k = 2:
 *   Acl = a − bk = 1 − 2 = −1 (Hurwitz)
 *
 * Under periodic sampling: control updates every T_s seconds.
 * Under event-triggering: control updates only when |e(t)| > σ|x(t)|.
 *
 * This example demonstrates:
 *   1. The trade-off between σ (threshold) and number of events
 *   2. Communication savings of ETC vs periodic sampling
 *   3. Effect of σ on convergence rate
 * ============================================================================ */

int main(void) {
    printf("=== Example 1: First-Order Stabilization with ETC ===\n\n");

    /* System parameters: A = [a], B = [b], K = [-k] */
    double A[] = {1.0};  /* Unstable open-loop */
    double B[] = {1.0};
    double K[] = {-2.0}; /* Stabilizing feedback: A+BK = -1 */

    double x0[] = {10.0}; /* Initial state */
    double T = 10.0;      /* Simulation horizon */
    double dt = 0.001;    /* Integration step */

    /* Compare three threshold values */
    double sigmas[] = {0.01, 0.1, 0.5};
    int n_cases = 3;

    printf("%-10s %-12s %-12s %-12s %-12s %s\n",
           "sigma", "events", "min IET", "avg IET", "final |x|", "comm ratio");
    printf("----------------------------------------------------------------------\n");

    for (int c = 0; c < n_cases; c++) {
        ETCSystem* sys = etc_system_create(A, B, K, 1, 1);
        etc_system_set_initial_state(sys, x0);
        etc_system_set_trigger(sys, ETC_TRIGGER_STATIC, sigmas[c], 0.0,
                                etc_trigger_static);

        etc_system_simulate(sys, T, dt);

        etc_history_compute_stats(&sys->history);
        double comm_ratio = etc_communication_ratio(
            sys->event_count, sys->t, dt);

        printf("%-10.3f %-12d %-12.6f %-12.6f %-12.6f %.3f\n",
               sigmas[c],
               sys->event_count,
               sys->history.min_iet,
               sys->history.avg_iet,
               fabs(sys->x.data[0]),
               comm_ratio);

        /* Periodic baseline: events if sampled every dt */
        double periodic_samples = T / dt;
        printf("  Periodic equivalent: %.0f samples, ETC savings: %.1f%%\n",
               periodic_samples,
               100.0 * (1.0 - comm_ratio));

        etc_system_free(sys);
    }

    printf("\nKey observations:\n");
    printf("  - Smaller sigma → more events (tighter control, faster convergence)\n");
    printf("  - Larger sigma → fewer events (communication savings)\n");
    printf("  - ETC achieves significant reduction vs periodic sampling\n");
    printf("  - All sigma < sigma_max guarantee asymptotic stability\n");

    return 0;
}
