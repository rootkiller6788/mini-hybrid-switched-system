/**
 * @file bench_hss_sim.c
 * @brief Performance benchmark for HSS simulation engine
 *
 * Benchmarks the RK4 integrator, event detection, and full
 * hybrid simulation loop for various system sizes.
 */
#include "hss_simulation.h"
#include "hss_core.h"
#include <stdio.h>
#include <time.h>

static double get_time_ms(void) {
    return (double)clock() / (double)(CLOCKS_PER_SEC) * 1000.0;
}

int main(void) {
    printf("=== HSS Simulation Benchmarks ===\n\n");

    double t_start, t_end;

    /* Benchmark 1: Bouncing ball 1M steps */
    printf("Bench 1: Bouncing ball (1M steps)...\n");
    HSS_BouncingBall bb = hss_bouncing_ball_init(9.81, 0.9, 100.0, 0.0);
    t_start = get_time_ms();
    for (int i = 0; i < 1000000; i++) {
        hss_bouncing_ball_step(&bb, 1e-5);
        if (bb.at_rest) break;
    }
    t_end = get_time_ms();
    printf("  Steps simulated, time = %.2f ms, bounces = %d\n",
           t_end - t_start, bb.bounce_count);

    /* Benchmark 2: RK4 integration */
    printf("\nBench 2: RK4 integration (10k steps, dim=10)...\n");
    HSS_System *sys = hss_system_create("bench", 1, 10, 0);
    double *A = calloc(100, sizeof(double));
    for (int i = 0; i < 10; i++) A[i * 10 + i] = -1.0 - i * 0.5;
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    double x0[10];
    for (int i = 0; i < 10; i++) x0[i] = 1.0;
    hss_set_initial_state(sys, 0, x0);
    free(A);

    t_start = get_time_ms();
    for (int i = 0; i < 10000; i++) {
        hss_step_continuous(sys, 1e-3);
    }
    t_end = get_time_ms();
    printf("  Time: %.2f ms (%.2f μs/step)\n",
           t_end - t_start, (t_end - t_start) / 10.0);

    hss_system_destroy(sys);

    /* Benchmark 3: Thermostat */
    printf("\nBench 3: Thermostat simulation (100k steps)...\n");
    HSS_Thermostat t = hss_thermostat_init(20, 18, 22, 0.8, 0.3, 30, 10);
    t_start = get_time_ms();
    for (int i = 0; i < 100000; i++) {
        hss_thermostat_step(&t, 1e-4);
    }
    t_end = get_time_ms();
    printf("  Time: %.2f ms\n", t_end - t_start);

    /* Benchmark 4: Full hybrid simulation */
    printf("\nBench 4: Full hybrid simulation (2-mode switched)...\n");
    sys = hss_system_create("bench2", 2, 2, 0);
    double A1[4] = {-0.1, 1, -1, -0.1};
    double A2[4] = {-2, 0, 0, -3};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A1, NULL, NULL);
    hss_mode_set_dynamics(sys, 1, HSS_CLASS_LINEAR, A2, NULL, NULL);

    HSS_Guard g;
    double ga[2] = {1, 0};
    g.a = ga; g.dim = 2; g.b = 5.0; g.is_active = true;
    hss_add_transition(sys, 0, 1, &g, 1, HSS_RESET_IDENTITY, NULL, NULL, "t");

    x0[0] = 1.0; x0[1] = 0.0;
    hss_set_initial_state(sys, 0, x0);

    HSS_SimConfig cfg = hss_sim_config_default();
    cfg.t_end = 10.0;
    cfg.dt = 0.001;
    cfg.max_steps = 20000;
    HSS_ExecutionTrace *trace = hss_trace_alloc(20000, 2);
    HSS_SimStats stats;

    t_start = get_time_ms();
    int ret = hss_simulate(sys, &cfg, trace, &stats);
    t_end = get_time_ms();

    printf("  Result: %d\n", ret);
    printf("  Time: %.2f ms\n", t_end - t_start);
    printf("  Steps: %d, Events: %d, Avg step: %.4f\n",
           stats.num_steps, stats.num_events, stats.avg_step);

    hss_trace_free(trace);
    hss_system_destroy(sys);

    printf("\n=== Benchmarks Complete ===\n");
    return 0;
}
