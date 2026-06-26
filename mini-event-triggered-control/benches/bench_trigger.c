#include "etc_core.h"
#include "etc_trigger.h"
#include "etc_dynamics.h"
#include <math.h>
#include <stdio.h>
#include <time.h>

/* ============================================================================
 * Benchmark: Trigger Evaluation Throughput
 *
 * Measures the computational cost of event-triggered control primitives:
 *  - Trigger function evaluation
 *  - Matrix exponential computation
 *  - Flow map evaluation
 *  - Simulation throughput
 * ============================================================================ */

static double wall_seconds(void) {
    return (double)clock() / (double)CLOCKS_PER_SEC;
}

int main(void) {
    printf("=== ETC Performance Benchmarks ===\n\n");

    double A[] = {0.0, 1.0, 0.0, 0.0};
    double B[] = {0.0, 1.0};
    double K[] = {-1.0, -2.0};
    double x0[] = {5.0, 0.0};

    const int N = 100000;
    double start, elapsed;

    /* ---- Trigger evaluation benchmark ---- */
    printf("--- Trigger Evaluation (%d iterations) ---\n", N);
    ETCVector x = etc_vector_create(2);
    ETCVector e = etc_vector_create(2);
    x.data[0] = 1.5; x.data[1] = 0.3;
    e.data[0] = 0.1; e.data[1] = 0.05;

    start = wall_seconds();
    volatile double sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += etc_trigger_static(&x, &e, 0.1, 0.0);
        sum += etc_trigger_quadratic(&x, &e, 0.1, 0.0);
        sum += etc_trigger_absolute(&x, &e, 0.0, 0.2);
    }
    elapsed = wall_seconds() - start;
    printf("  3 × %d trigger evaluations: %.4f s (%.2f M evals/s)\n",
           N, elapsed, (3.0 * N) / (elapsed * 1e6 + 1e-12));
    printf("  (checksum: %.2f)\n", sum);

    etc_vector_free(&x); etc_vector_free(&e);

    /* ---- Matrix exponential benchmark ---- */
    printf("\n--- Matrix Exponential (%d iterations) ---\n", N / 100);
    ETCMatrix M = etc_matrix_create(2, 2);
    M.data[0] = 0.0; M.data[1] = 1.0;
    M.data[2] = -1.0; M.data[3] = -2.0;
    ETCMatrix expM = etc_matrix_create(2, 2);

    start = wall_seconds();
    for (int i = 0; i < N / 100; i++) {
        etc_matrix_exponential(&M, 0.1 * (double)(i % 10 + 1), 2, &expM);
    }
    elapsed = wall_seconds() - start;
    printf("  %d 2×2 exp(A t): %.4f s (%.2f K exp/s)\n",
           N / 100, elapsed, (N / 100) / (elapsed * 1e3 + 1e-12));

    etc_matrix_free(&M); etc_matrix_free(&expM);

    /* ---- Simulation throughput benchmark ---- */
    printf("\n--- Full Simulation Throughput ---\n");
    ETCSystem* sys = etc_system_create(A, B, K, 2, 1);
    etc_system_set_initial_state(sys, x0);
    etc_system_set_trigger(sys, ETC_TRIGGER_STATIC, 0.05, 0.0,
                            etc_trigger_static);

    start = wall_seconds();
    etc_system_simulate(sys, 10.0, 0.001);
    elapsed = wall_seconds() - start;

    printf("  10s simulation (dt=0.001, %d steps, %d events): %.4f s\n",
           (int)(10.0 / 0.001), sys->event_count, elapsed);
    printf("  Real-time factor: %.2f×\n", 10.0 / (elapsed + 1e-12));
    printf("  Avg step time: %.3f μs\n",
           elapsed * 1e6 / (10.0 / 0.001));

    etc_system_free(sys);

    /* ---- IET computation benchmark ---- */
    printf("\n--- IET Lower Bound Computation ---\n");
    sys = etc_system_create(A, B, K, 2, 1);
    start = wall_seconds();
    for (int i = 0; i < N; i++) {
        volatile double bnd = etc_inter_event_time_bound(sys);
        (void)bnd;
    }
    elapsed = wall_seconds() - start;
    printf("  %d IET bound computations: %.4f s\n", N, elapsed);
    etc_system_free(sys);

    printf("\n=== Benchmarks Complete ===\n");
    return 0;
}
