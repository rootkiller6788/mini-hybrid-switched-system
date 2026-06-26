/**
 * @file bouncing_ball.c
 * @brief Bouncing ball example — canonical Zeno hybrid automaton
 *
 * Demonstrates execution of the bouncing ball, which exhibits
 * Zeno behavior: infinitely many bounces in finite time.
 */

#include "../include/hybrid_automaton.h"
#include "../include/hybrid_execution.h"
#include "../include/hybrid_examples.h"
#include "../include/hybrid_simulation.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("=== Bouncing Ball Hybrid Automaton ===\n\n");

    /* Create bouncing ball: h0=10m, v0=0, c=0.8, g=9.81 */
    double h0 = 10.0, v0 = 0.0, c = 0.8, g = 9.81;
    HybridAutomaton *ha = example_bouncing_ball(h0, v0, c, g);

    if (!ha) {
        printf("Failed to create automaton\n");
        return 1;
    }

    printf("System: %s\n", ha->name);
    printf("Variables: height (h), velocity (v)\n");
    printf("Initial: h0=%.2f m, v0=%.2f m/s\n", h0, v0);
    printf("Restitution: c=%.2f\n", c);
    printf("Gravity: g=%.2f m/s²\n\n", g);

    /* Print automaton structure */
    hybrid_automaton_print(ha);

    /* Simulate */
    HybridSimConfig config = HYBRID_SIM_CONFIG_DEFAULT;
    config.dt = 0.001;
    config.t_max = 5.0;
    config.max_transitions = 100;
    config.stop_on_zeno = true;
    config.verbosity = 0;

    printf("\n--- Simulation ---\n");
    HybridExecution *exec = hybrid_simulate(ha, &config);

    if (exec) {
        printf("Simulation completed:\n");
        printf("  Total time: %.4f s\n", exec->total_time);
        printf("  Flow segments: %d\n", exec->num_segments);
        printf("  Discrete jumps (bounces): %d\n", exec->num_jumps);
        printf("  Zeno detected: %s\n", exec->is_zeno ? "YES" : "NO");

        if (exec->is_zeno) {
            double T_inf = hybrid_execution_zeno_time(exec);
            printf("  Estimated Zeno time: %.4f s\n", T_inf);

            /* Theoretical: T∞ = v0/g + 2v0/(g(1-c)) for initial v0=0
               After first fall: v1 = sqrt(2*g*h0)
               T∞ = sqrt(2h0/g) * (1 + 2c/(1-c)) */
            double v1 = sqrt(2.0 * g * h0);
            double T_theory = v1 / g + 2.0 * v1 / (g * (1.0 - c));
            printf("  Theoretical Zeno time: %.4f s\n", T_theory);
        }

        /* Print first 5 bounces */
        int show = exec->num_jumps < 5 ? exec->num_jumps : 5;
        printf("\n  First %d bounces:\n", show);
        for (int i = 0; i < show; i++) {
            printf("    t=%.4f: h=%.4f, v=%.4f → v'=%.4f\n",
                   exec->jump_points[i].t_jump,
                   exec->jump_points[i].x_pre[0],
                   exec->jump_points[i].x_pre[1],
                   exec->jump_points[i].x_post[1]);
        }

        hybrid_execution_destroy(exec);
    } else {
        printf("Simulation failed\n");
    }

    hybrid_automaton_destroy(ha);
    return 0;
}
