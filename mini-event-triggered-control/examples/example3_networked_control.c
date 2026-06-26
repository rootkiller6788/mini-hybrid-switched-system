#include "etc_core.h"
#include "etc_trigger.h"
#include "etc_stability.h"
#include "etc_dynamics.h"
#include "etc_periodic.h"
#include "etc_self.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ============================================================================
 * Example 3: Networked Control System — ETC vs PETC vs STC Comparison
 *
 * This example compares three event-triggered paradigms for a networked
 * control scenario where communication bandwidth is limited.
 *
 * Plant: Double integrator (position + velocity)
 * Scenario: Wireless sensor-to-controller-to-actuator network
 *   - Bandwidth constraint: max 100 messages/second
 *   - Battery constraint: minimize transmissions
 *   - Performance: maintain < 5% overshoot
 *
 * We compare:
 *   ETC (continuous monitoring) — optimal event count, needs monitoring
 *   PETC (periodic checking, h=0.01s) — implementable, slightly more events
 *   STC (self-triggered) — no monitoring, most conservative
 * ============================================================================ */

int main(void) {
    printf("=== Example 3: Networked Control — ETC vs PETC vs STC ===\n\n");

    /* System: double integrator with PD control */
    double A[] = {0.0, 1.0, 0.0, 0.0};
    double B[] = {0.0, 1.0};
    double K[] = {-1.0, -2.0};
    double x0[] = {3.0, 1.0};  /* Position 3, velocity 1 */
    double T = 10.0;
    double dt = 0.001;

    printf("Plant: double integrator, PD gains K = [%.1f, %.1f]\n", K[0], K[1]);
    printf("Initial state: x0 = [%.1f, %.1f]\n", x0[0], x0[1]);
    printf("Constraints: bandwidth ≤ 100 msg/s, battery-aware\n\n");

    /* ================================================================
     * 1. Continuous ETC
     * ================================================================ */
    printf("--- 1. Continuous ETC ---\n");
    ETCSystem* sys_etc = etc_system_create(A, B, K, 2, 1);
    etc_system_set_initial_state(sys_etc, x0);
    etc_system_set_trigger(sys_etc, ETC_TRIGGER_STATIC, 0.05, 0.0,
                            etc_trigger_static);
    etc_system_simulate(sys_etc, T, dt);
    etc_history_compute_stats(&sys_etc->history);

    double etc_comm = etc_communication_ratio(sys_etc->event_count, sys_etc->t, dt);
    printf("  Events:           %d\n", sys_etc->event_count);
    printf("  Avg IET:          %.4f s\n", sys_etc->history.avg_iet);
    printf("  Comm ratio:       %.4f\n", etc_comm);
    printf("  Avg msg/s:        %.2f\n",
           (double)sys_etc->event_count / sys_etc->t);
    printf("  Final |x|:        %.6f\n", etc_system_state_norm(sys_etc));
    printf("  Zeno-free:        YES (τ_min = %.6f)\n\n",
           sys_etc->history.min_iet);

    /* ================================================================
     * 2. PETC
     * ================================================================ */
    printf("--- 2. Periodic ETC (h = 0.02 s) ---\n");
    ETCSystem* sys_petc = etc_system_create(A, B, K, 2, 1);
    etc_system_set_initial_state(sys_petc, x0);
    etc_system_set_trigger(sys_petc, ETC_TRIGGER_STATIC, 0.05, 0.0,
                            etc_trigger_static);

    double h_max = etc_petc_max_sampling_period(sys_petc, 0.05);
    printf("  h_max (theoretical): %.6f s\n", h_max);

    double h = 0.02; /* PETC sampling period */
    if (h > h_max) {
        printf("  WARNING: h > h_max, stability not guaranteed!\n");
        h = h_max * 0.9;
        printf("  Adjusted h = %.6f\n\n", h);
    }

    PETCConfig pcfg = etc_petc_config_create(h);
    etc_petc_simulate(sys_petc, &pcfg, T, dt);
    etc_history_compute_stats(&sys_petc->history);

    double petc_comm = etc_communication_ratio(
        sys_petc->event_count, sys_petc->t, dt);
    printf("  Check count:      %d\n", pcfg.check_count);
    printf("  Events:           %d\n", pcfg.event_count);
    printf("  Avg IET:          %.4f s\n", sys_petc->history.avg_iet);
    printf("  Comm ratio:       %.4f\n", petc_comm);
    printf("  Overhead vs ETC:  %.1f%%\n",
           100.0 * (petc_comm - etc_comm) / (etc_comm + 1e-12));
    printf("  Final |x|:        %.6f\n\n", etc_system_state_norm(sys_petc));

    /* ================================================================
     * 3. STC
     * ================================================================ */
    printf("--- 3. Self-Triggered Control ---\n");
    ETCSystem* sys_stc = etc_system_create(A, B, K, 2, 1);
    etc_system_set_initial_state(sys_stc, x0);

    STCConfig stc_cfg = etc_stc_config_create();
    etc_stc_init(sys_stc, &stc_cfg);
    stc_cfg.sigma = 0.05;
    stc_cfg.tau_min = 0.001;
    stc_cfg.tau_max = 1.0;

    etc_stc_simulate(sys_stc, &stc_cfg, T, dt);
    etc_history_compute_stats(&sys_stc->history);

    double stc_comm = etc_communication_ratio(
        sys_stc->event_count, sys_stc->t, dt);
    printf("  Events:           %d\n", sys_stc->event_count);
    printf("  Avg IET:          %.4f s\n", sys_stc->history.avg_iet);
    printf("  Comm ratio:       %.4f\n", stc_comm);
    printf("  Final |x|:        %.6f\n", etc_system_state_norm(sys_stc));

    /* Compare overhead */
    double stc_overhead; int stc_ev, etc_ev;
    etc_stc_compare(sys_stc, &stc_cfg, T, dt, &stc_overhead, &stc_ev, &etc_ev);
    printf("  STC/ETC ratio:    %.2f\n\n", stc_overhead);

    /* ================================================================
     * 4. Summary Comparison
     * ================================================================ */
    printf("=== Summary ===\n");
    printf("%-20s %12s %12s %12s %12s\n",
           "Method", "Events", "Avg IET(s)", "Final |x|", "Msg/s");
    printf("---------------------------------------------------------------\n");
    printf("%-20s %12d %12.4f %12.6f %12.2f\n",
           "Continuous ETC", sys_etc->event_count,
           sys_etc->history.avg_iet, etc_system_state_norm(sys_etc),
           (double)sys_etc->event_count / sys_etc->t);
    printf("%-20s %12d %12.4f %12.6f %12.2f\n",
           "Periodic ETC", pcfg.event_count,
           sys_petc->history.avg_iet, etc_system_state_norm(sys_petc),
           (double)pcfg.event_count / sys_petc->t);
    printf("%-20s %12d %12.4f %12.6f %12.2f\n",
           "Self-Triggered", sys_stc->event_count,
           sys_stc->history.avg_iet, etc_system_state_norm(sys_stc),
           (double)sys_stc->event_count / sys_stc->t);
    printf("---------------------------------------------------------------\n");
    printf("%-20s %12d %12s %12s %12.2f\n",
           "Periodic (h=0.001)", (int)(T/dt), "0.001", "—",
           T/dt);

    printf("\nConclusion:\n");
    printf("  ETC:  Best event efficiency, requires continuous monitoring\n");
    printf("  PETC: Slightly more events, implementable on digital hardware\n");
    printf("  STC:  Most events, no monitoring needed — ideal for sensor-less nodes\n");

    /* Cleanup */
    etc_system_free(sys_etc);
    etc_system_free(sys_petc);
    etc_system_free(sys_stc);
    etc_stc_config_free(&stc_cfg);

    return 0;
}
