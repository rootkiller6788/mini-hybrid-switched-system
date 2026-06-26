#include "etc_core.h"
#include "etc_trigger.h"
#include "etc_stability.h"
#include "etc_dynamics.h"
#include "etc_periodic.h"
#include "etc_self.h"
#include <math.h>
#include <stdio.h>

/* ============================================================================
 * Demo: Event-Triggered Control Overview
 *
 * Interactive demonstration of key ETC concepts:
 *   1. Trigger function visualization
 *   2. Event pattern analysis
 *   3. Parameter sensitivity (σ, h)
 * ============================================================================ */

static void print_separator(const char* title) {
    printf("\n========================================\n");
    printf("  %s\n", title);
    printf("========================================\n\n");
}

static void demonstrate_trigger_mechanics(void) {
    print_separator("1. Trigger Function Mechanics");

    ETCVector x = etc_vector_create(2);
    ETCVector e = etc_vector_create(2);
    x.data[0] = 2.0; x.data[1] = 0.0;

    printf("State x = [%.1f, %.1f], threshold σ = 0.1\n", x.data[0], x.data[1]);
    printf("\n%-12s %-12s %-12s %s\n",
           "|e|", "|x|", "Γ(e)", "Trigger?");
    printf("------------------------------------------\n");

    for (int i = 0; i <= 10; i++) {
        double error = 0.2 * (double)i / 10.0;
        e.data[0] = error; e.data[1] = 0.0;
        double gamma = etc_trigger_static(&x, &e, 0.1, 0.0);
        printf("%-12.6f %-12.6f %-12.6f %s\n",
               error, etc_vector_norm(&x), gamma,
               gamma >= 0.0 ? "YES ***" : "no");
    }

    etc_vector_free(&x); etc_vector_free(&e);
}

static void demonstrate_trigger_types(void) {
    print_separator("2. Trigger Type Comparison");

    ETCVector x = etc_vector_create(2);
    ETCVector e = etc_vector_create(2);
    x.data[0] = 1.0; x.data[1] = 0.0;
    e.data[0] = 0.1; e.data[1] = 0.1;

    double e_norm = etc_vector_norm(&e);
    double x_norm = etc_vector_norm(&x);

    printf("|x| = %.4f, |e| = %.4f, σ = 0.1, ε = 0.05\n\n", x_norm, e_norm);
    printf("%-18s %-12s %-12s %s\n", "Trigger Type", "Γ value", "Threshold", "Fires?");
    printf("-------------------------------------------------------\n");

    double gs = etc_trigger_static(&x, &e, 0.1, 0.0);
    printf("%-18s %-12.6f %-12s %s\n", "Static", gs, "σ|x| = 0.1",
           gs >= 0.0 ? "YES" : "no");

    double gq = etc_trigger_quadratic(&x, &e, 0.1, 0.0);
    printf("%-18s %-12.6f %-12s %s\n", "Quadratic", gq, "σ|x|² = 0.1",
           gq >= 0.0 ? "YES" : "no");

    double ga = etc_trigger_absolute(&x, &e, 0.0, 0.3);
    printf("%-18s %-12.6f %-12s %s\n", "Absolute", ga, "ε = 0.3",
           ga >= 0.0 ? "YES" : "no");

    double gm = etc_trigger_mixed(&x, &e, 0.1, 0.02);
    printf("%-18s %-12.6f %-12s %s\n", "Mixed", gm, "σ|x|²+ε = 0.12",
           gm >= 0.0 ? "YES" : "no");

    etc_vector_free(&x); etc_vector_free(&e);
}

static void demonstrate_parameter_sweep(void) {
    print_separator("3. Parameter Sensitivity (σ sweep)");

    double A[] = {0.0, 1.0, 0.0, 0.0};
    double B[] = {0.0, 1.0};
    double K[] = {-1.0, -2.0};
    double x0[] = {5.0, 0.0};
    double T = 8.0, dt = 0.001;

    printf("Double integrator, PD control, T=%.0fs\n\n", T);
    printf("%-10s %-10s %-12s %-12s %-12s\n",
           "σ", "Events", "Min IET", "Avg IET", "Final |x|");
    printf("----------------------------------------------------------\n");

    double sigmas[] = {0.01, 0.02, 0.05, 0.1, 0.2, 0.5};
    for (int i = 0; i < 6; i++) {
        ETCSystem* sys = etc_system_create(A, B, K, 2, 1);
        etc_system_set_initial_state(sys, x0);
        etc_system_set_trigger(sys, ETC_TRIGGER_STATIC, sigmas[i], 0.0,
                                etc_trigger_static);
        etc_system_simulate(sys, T, dt);
        etc_history_compute_stats(&sys->history);
        printf("%-10.3f %-10d %-12.6f %-12.6f %-12.6f\n",
               sigmas[i], sys->event_count,
               sys->history.min_iet, sys->history.avg_iet,
               etc_system_state_norm(sys));
        etc_system_free(sys);
    }
}

static void demonstrate_etc_vs_time_triggered(void) {
    print_separator("4. ETC vs Time-Triggered Communication");

    double A[] = {0.0, 1.0, 0.0, 0.0};
    double B[] = {0.0, 1.0};
    double K[] = {-1.0, -2.0};
    double x0[] = {5.0, 0.0};
    double T = 10.0, dt = 0.001;

    printf("Time-Triggered (periodic) vs Event-Triggered\n\n");
    printf("%-18s %12s %12s %12s\n",
           "Method", "Updates", "Avg IET(s)", "Final |x|");
    printf("----------------------------------------------------------\n");

    /* Time-triggered with different periods */
    double periods[] = {0.01, 0.05, 0.1, 0.5};
    for (int i = 0; i < 4; i++) {
        int updates = (int)(T / periods[i]) + 1;
        printf("%-18s %12d %12.4f %12s\n",
               "Periodic", updates, periods[i], "—");
    }

    /* ETC with different σ */
    double sigmas[] = {0.01, 0.05, 0.1, 0.5};
    for (int i = 0; i < 4; i++) {
        ETCSystem* sys = etc_system_create(A, B, K, 2, 1);
        etc_system_set_initial_state(sys, x0);
        etc_system_set_trigger(sys, ETC_TRIGGER_STATIC, sigmas[i], 0.0,
                                etc_trigger_static);
        etc_system_simulate(sys, T, dt);
        etc_history_compute_stats(&sys->history);
        printf("%-18s %12d %12.4f %12.6f\n",
               "Event-Triggered", sys->event_count,
               sys->history.avg_iet, etc_system_state_norm(sys));
        etc_system_free(sys);
    }

    printf("\nKey insight: ETC adapts update rate to state evolution.\n");
    printf("When the system is near equilibrium, few updates are needed.\n");
    printf("When the system is rapidly changing, updates occur more frequently.\n");
}

int main(void) {
    printf("╔══════════════════════════════════════════╗\n");
    printf("║  Event-Triggered Control — Demo Tour    ║\n");
    printf("║  Based on Tabuada (2007), Heemels (2012)║\n");
    printf("╚══════════════════════════════════════════╝\n");

    demonstrate_trigger_mechanics();
    demonstrate_trigger_types();
    demonstrate_parameter_sweep();
    demonstrate_etc_vs_time_triggered();

    print_separator("End of Demo");
    printf("See examples/ and tests/ for more detailed usage.\n\n");
    return 0;
}
