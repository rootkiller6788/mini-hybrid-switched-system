#include "etc_core.h"
#include "etc_trigger.h"
#include "etc_stability.h"
#include "etc_dynamics.h"
#include <math.h>
#include <stdio.h>

/* ============================================================================
 * Example 2: Double Integrator with Event-Triggered Control
 *
 * Plant:  ẋ₁ = x₂,  ẋ₂ = u   (double integrator, open-loop unstable)
 * Controller: u = −k₁ x₁ − k₂ x₂  (PD-like state feedback)
 *
 * With k₁ = 1, k₂ = 2:
 *   Acl = [0 1; −1 −2], eigenvalues at −1 (double), critically damped.
 *
 * This example demonstrates:
 *   1. State-dependent inter-event times for 2D system
 *   2. Lyapunov function computation
 *   3. ISS condition verification
 *   4. Event-triggering pattern visualization (text-based)
 *   5. Trigger margin evolution over time
 * ============================================================================ */

int main(void) {
    printf("=== Example 2: Double Integrator with ETC ===\n\n");

    /* System matrices */
    double A[] = {0.0, 1.0, 0.0, 0.0};
    double B[] = {0.0, 1.0};
    double K[] = {-1.0, -2.0};  /* PD gains */
    double x0[] = {5.0, 0.0};   /* Initial: position=5, velocity=0 */

    ETCSystem* sys = etc_system_create(A, B, K, 2, 1);
    etc_system_set_initial_state(sys, x0);
    etc_system_set_trigger(sys, ETC_TRIGGER_STATIC, 0.05, 0.0,
                            etc_trigger_static);

    /* Compute the Control Lyapunov Function */
    printf("Computing CLF (Acl^T P + P Acl = -I)...\n");
    if (etc_compute_clf(sys)) {
        printf("  P = [%.4f %.4f; %.4f %.4f]\n",
               sys->V.P.data[0], sys->V.P.data[1],
               sys->V.P.data[2], sys->V.P.data[3]);
        printf("  λ_min(P) = %.4f, λ_max(P) = %.4f\n",
               sys->V.lambda_min_P, sys->V.lambda_max_P);
        printf("  P > 0: YES (Hurwitz closed-loop)\n");
    } else {
        printf("  CLF computation failed (check Acl eigenvalues)\n");
    }

    /* ISS Lyapunov verification */
    double a1, a2, a3, g;
    if (etc_verify_iss_lyapunov(sys, &a1, &a2, &a3, &g)) {
        printf("  ISS Lyapunov: VERIFIED\n");
        printf("    α₁(s) = %.4f·s², α₂(s) = %.4f·s²\n", a1, a2);
        printf("    α₃(s) = %.4f·s², γ(s) = %.4f·s²\n", a3, g);
        printf("    Dissipation margin: α₃(|x|) > γ(|e|) for |e| < %.4f|x|\n",
               sqrt(a3 / (g + 1e-12)));
    } else {
        printf("  ISS Lyapunov: NOT VERIFIED (||PBK|| too large)\n");
    }

    /* Compute σ_max and verify our σ */
    ETCMatrix BK = etc_matrix_create(2, 2);
    etc_matrix_mul(&sys->B, &sys->K, &BK);
    double sigma_max = etc_compute_sigma_max(&sys->V.P, &BK);
    etc_matrix_free(&BK);
    printf("  σ_max = %.6f (current σ=%.3f, stabilizing: %s)\n",
           sigma_max, sys->sigma,
           etc_is_sigma_stabilizing(&sys->V.P, &BK, sys->sigma) ? "YES" : "NO");

    /* Theoretical IET bound */
    double tau_bound = etc_inter_event_time_bound(sys);
    printf("  Theoretical τ_min ≥ %.6f s\n\n", tau_bound);

    /* Simulate */
    double T = 15.0, dt = 0.001;
    printf("Simulating T = %.0f seconds...\n", T);
    etc_system_simulate(sys, T, dt);

    /* Results */
    etc_history_compute_stats(&sys->history);
    printf("\n--- Results ---\n");
    printf("Events triggered:       %d\n", sys->event_count);
    printf("Min inter-event time:   %.6f s\n", sys->history.min_iet);
    printf("Max inter-event time:   %.6f s\n", sys->history.max_iet);
    printf("Avg inter-event time:   %.6f s\n", sys->history.avg_iet);
    printf("Final state:            [%.6f, %.6f]\n",
           sys->x.data[0], sys->x.data[1]);
    printf("Final |x|:              %.6f\n", etc_system_state_norm(sys));
    printf("Communication ratio:    %.3f\n",
           etc_communication_ratio(sys->event_count, sys->t, dt));
    printf("Zeno-free:              %s\n",
           etc_check_zeno_free(sys) ? "YES" : "NO");

    /* Print event timeline (first 15 events) */
    printf("\n--- Event Timeline (first 15) ---\n");
    printf("%-6s %-10s %-10s %-12s %-12s\n",
           "k", "t_k", "Δt_k", "|x(t_k)|", "Γ(t_k)");
    for (int i = 0; i < sys->history.n_events && i < 15; i++) {
        ETCEvent* ev = &sys->history.events[i];
        printf("%-6d %-10.4f %-10.6f %-12.6f %-12.6f\n",
               ev->event_index, ev->time, ev->inter_event_time,
               etc_vector_norm(&ev->state), ev->trigger_value);
    }

    /* Practical stability check */
    double radius;
    if (etc_check_practical_stability(sys, &radius)) {
        printf("\nPractically stable: YES (radius = %.6f)\n", radius);
    } else {
        printf("\nPractically stable: NO (state outside bound)\n");
    }

    etc_system_free(sys);
    return 0;
}
