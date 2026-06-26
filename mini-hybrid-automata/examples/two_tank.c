/**
 * @file two_tank.c
 * @brief Two-tank system example — controlled water flow
 *
 * Demonstrates the two-tank hybrid automaton where a valve
 * controls inter-tank water flow. Modes: OPEN and CLOSED.
 */

#include "../include/hybrid_automaton.h"
#include "../include/hybrid_execution.h"
#include "../include/hybrid_examples.h"
#include "../include/hybrid_simulation.h"
#include <stdio.h>

int main(void)
{
    printf("=== Two-Tank Hybrid Automaton ===\n\n");

    double K_inter = 0.3, K_out = 0.05;
    double h1_init = 8.0, h2_init = 2.0;

    HybridAutomaton *ha = example_two_tank(K_inter, K_out, h1_init, h2_init);
    if (!ha) { printf("Failed\n"); return 1; }

    printf("Inter-tank coefficient: %.2f\n", K_inter);
    printf("Output coefficient: %.3f\n", K_out);
    printf("Initial: h1=%.1f, h2=%.1f (CLOSED)\n\n", h1_init, h2_init);

    hybrid_automaton_print(ha);

    /* Simulate with OPEN valve (we start CLOSED, so tank 2 drains) */
    HybridSimConfig config = HYBRID_SIM_CONFIG_DEFAULT;
    config.dt = 0.01;
    config.t_max = 30.0;
    config.max_transitions = 50;

    HybridExecution *exec = hybrid_simulate(ha, &config);
    if (exec) {
        printf("\n--- Simulation (starts CLOSED, drains tank 2) ---\n");
        printf("Total time: %.1f s\n", exec->total_time);
        printf("Flow segments: %d\n", exec->num_segments);

        if (exec->num_segments > 0) {
            /* Show final state of first segment */
            int last = exec->num_segments - 1;
            printf("Final state: h1=%.3f, h2=%.3f\n",
                   exec->flow_segments[last].x_end[0],
                   exec->flow_segments[last].x_end[1]);
        }

        hybrid_execution_destroy(exec);
    }

    hybrid_automaton_destroy(ha);
    return 0;
}
