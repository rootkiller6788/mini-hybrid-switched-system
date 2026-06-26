#include "dta_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* bench_core.c - Performance benchmark for core operations */

static double get_time_ms(clock_t start, clock_t end) {
    return 1000.0 * (double)(end - start) / (double)CLOCKS_PER_SEC;
}

int main(void) {
    printf("=== Dwell-Time Analysis Benchmarks ===\n\n");
    int n = 4;
    int n_runs = 1000;

    double A[16] = {-1,1,0,0, -1,-2,1,0, 0,-1,-3,1, 0,0,-1,-4};
    double re[4], im[4];

    /* Benchmark eigenvalue computation */
    clock_t t0 = clock();
    int i;
    for (i = 0; i < n_runs; i++)
        dta_eigenvalues(A, n, re, im, 200);
    clock_t t1 = clock();
    printf("Eigenvalues (%dx%d, %d runs): %.2f ms total, %.4f ms/op\n",
           n, n, n_runs, get_time_ms(t0, t1), get_time_ms(t0, t1)/n_runs);

    /* Benchmark matrix exponential */
    double exp_result[16];
    t0 = clock();
    for (i = 0; i < n_runs; i++)
        dta_matrix_exp(A, n, 0.1, exp_result);
    t1 = clock();
    printf("Matrix exp  (%dx%d, %d runs): %.2f ms total, %.4f ms/op\n",
           n, n, n_runs, get_time_ms(t0, t1), get_time_ms(t0, t1)/n_runs);

    /* Benchmark Lyapunov solver */
    double Q[16] = {2,0,0,0,0,2,0,0,0,0,2,0,0,0,0,2};
    double P[16];
    t0 = clock();
    for (i = 0; i < n_runs; i++)
        dta_solve_lyapunov(A, n, Q, P);
    t1 = clock();
    printf("Lyapunov    (%dx%d, %d runs): %.2f ms total, %.4f ms/op\n",
           n, n, n_runs, get_time_ms(t0, t1), get_time_ms(t0, t1)/n_runs);

    /* Benchmark Hurwitz check */
    t0 = clock();
    for (i = 0; i < n_runs; i++)
        dta_is_hurwitz(A, n, 1e-9);
    t1 = clock();
    printf("Hurwitz     (%dx%d, %d runs): %.2f ms total, %.4f ms/op\n",
           n, n, n_runs, get_time_ms(t0, t1), get_time_ms(t0, t1)/n_runs);

    printf("\n=== Benchmarks Complete ===\n");
    return 0;
}
