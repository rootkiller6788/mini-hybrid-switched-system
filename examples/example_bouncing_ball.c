/**
 * @file example_bouncing_ball.c
 * @brief L6 KP1: Canonical Bouncing Ball — Hybrid Automaton Example
 *
 * Demonstrates the classic hybrid system: a ball bouncing under gravity
 * with energy loss at each impact. This is the "hello world" of hybrid
 * systems theory, exhibiting Zeno behavior (infinitely many bounces
 * in finite time).
 *
 * Model:
 *   Flow:    ḣ = v, v̇ = -g           (free fall)
 *   Guard:   h ≤ 0 AND v < 0         (ground impact)
 *   Reset:   v⁺ = -c · v⁻            (restitution)
 *
 * Zeno time: τ_z = √(2h₀/g) · (1+c)/(1-c)
 * Total bounces: infinite (theoretically), energy decays geometrically.
 */

#include "hss_simulation.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Bouncing Ball — Canonical Hybrid Automaton (L6) ===\n\n");

    /* Initialize with 10m drop, 80% restitution */
    HSS_BouncingBall bb = hss_bouncing_ball_init(9.81, 0.8, 10.0, 0.0);

    printf("Parameters:\n");
    printf("  Gravity g = %.2f m/s²\n", bb.gravity);
    printf("  Restitution c = %.2f\n", bb.restitution);
    printf("  Initial height = %.2f m\n", bb.height);
    printf("  Initial energy = %.4f J/kg\n\n", bb.energy);

    /* Compute theoretical Zeno time */
    double t_first = sqrt(2.0 * 10.0 / 9.81);
    double zeno_time = t_first * (1.0 + bb.restitution)
                       / (1.0 - bb.restitution);
    printf("  Theoretical Zeno time = %.4f s\n", zeno_time);
    printf("  First impact time = %.4f s\n\n", t_first);

    /* Simulate step by step, printing bounce events */
    printf("Simulation trace:\n");
    printf("  %8s  %10s  %10s  %10s\n", "Time(s)", "Height(m)", "Vel(m/s)", "Bounces");

    double dt = 0.01;
    double max_time = 20.0;
    double print_interval = 0.5;
    double next_print = 0.0;

    while (bb.time < max_time && !bb.at_rest) {
        int prev_bounces = bb.bounce_count;
        hss_bouncing_ball_step(&bb, dt);

        if (bb.bounce_count > prev_bounces || bb.time >= next_print) {
            printf("  %8.4f  %10.6f  %10.6f  %10d\n",
                   bb.time, bb.height, bb.velocity, bb.bounce_count);
            next_print += print_interval;
        }

        if (bb.bounce_count > 50) {
            printf("  ... (Zeno accumulation: stopped at 50 bounces)\n");
            break;
        }
    }

    printf("\nResults:\n");
    printf("  Total bounces: %d\n", bb.bounce_count);
    printf("  Final state: h=%.6f m, v=%.6f m/s, at_rest=%s\n",
           bb.height, bb.velocity, bb.at_rest ? "yes" : "no");
    printf("  Total simulation time: %.4f s\n", bb.time);

    /* Verify energy decay */
    double initial_energy = 9.81 * 10.0;  /* mgh₀ */
    double final_energy = 9.81 * bb.height
                         + 0.5 * bb.velocity * bb.velocity;
    printf("  Energy retention: %.2f%%\n",
           100.0 * final_energy / initial_energy);

    return 0;
}
