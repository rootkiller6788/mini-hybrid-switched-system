/**
 * @file thermostat.c
 * @brief Thermostat example — mode-switching temperature control
 *
 * Demonstrates the classic thermostat hybrid automaton where
 * heating toggles ON/OFF to maintain temperature in a target range.
 */

#include "../include/hybrid_automaton.h"
#include "../include/hybrid_execution.h"
#include "../include/hybrid_examples.h"
#include "../include/hybrid_simulation.h"
#include <stdio.h>

int main(void)
{
    printf("=== Thermostat Hybrid Automaton ===\n\n");

    double T_lo = 18.0, T_hi = 22.0, T_env = 15.0;
    double K_h = 0.05, K_on = 1.5, T_init = 20.0;

    HybridAutomaton *ha = example_thermostat(T_lo, T_hi, T_env,
                                               K_h, K_on, T_init, true);
    if (!ha) { printf("Failed\n"); return 1; }

    printf("Target range: [%.1f, %.1f]°C\n", T_lo, T_hi);
    printf("Environment: %.1f°C\n", T_env);
    printf("Heating rate: %.1f°C/s\n", K_on);
    printf("Initial: %.1f°C (ON)\n\n", T_init);

    hybrid_automaton_print(ha);

    HybridSimConfig config = HYBRID_SIM_CONFIG_DEFAULT;
    config.dt = 0.01;
    config.t_max = 60.0;
    config.max_transitions = 200;

    HybridExecution *exec = hybrid_simulate(ha, &config);
    if (exec) {
        printf("\n--- Simulation Results ---\n");
        printf("Total time: %.1f s\n", exec->total_time);
        printf("Mode switches: %d\n", exec->num_jumps);

        /* Show switching sequence */
        printf("\nSwitching schedule (first 10):\n");
        int show = exec->num_jumps < 10 ? exec->num_jumps : 10;
        for (int i = 0; i < show; i++) {
            const char *src = ha->modes[exec->jump_points[i].src_mode].name;
            const char *tgt = ha->modes[exec->jump_points[i].tgt_mode].name;
            printf("  t=%.2f: %s → %s, T=%.3f\n",
                   exec->jump_points[i].t_jump, src, tgt,
                   exec->jump_points[i].x_pre[0]);
        }

        hybrid_execution_destroy(exec);
    }

    hybrid_automaton_destroy(ha);
    return 0;
}
