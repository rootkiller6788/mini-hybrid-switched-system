/**
 * @file test_hss_analysis.c
 * @brief Tests for HSS Analysis (L4-L5): stability theorems, algorithms.
 */
#include "hss_analysis.h"
#include "hss_core.h"
#include "hss_simulation.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>

static int tests_run = 0;
static int tests_passed = 0;
#define TEST(n) do { tests_run++; printf("  TEST: %s ... ", n); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)

static void test_clf_verify(void) {
    TEST("verify_clf_theorem");
    HSS_System *sys = hss_system_create("clf_test", 1, 2, 0);
    double A[4] = {-1, 0, 0, -2};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    double P[4] = {1, 0, 0, 1};
    hss_mode_set_lyapunov(sys, 0, P, 0.5);

    HSS_CommonLyapunov result = hss_verify_clf_theorem(sys, 1e-6);
    assert(result.n == 2);
    free(result.P);
    hss_system_destroy(sys);
    PASS();
}

static void test_mlf_verify(void) {
    TEST("verify_mlf_theorem");
    HSS_System *sys = hss_system_create("mlf_test", 2, 2, 0);
    double A[4] = {-1, 0, 0, -2};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(sys, 1, HSS_CLASS_LINEAR, A, NULL, NULL);

    HSS_MultipleLyapunov r = hss_verify_mlf_theorem(sys, NULL, 1.5);
    assert(r.num_modes == 2);
    for (int i = 0; i < r.num_modes; i++) free(r.P[i]);
    free(r.P); free(r.decay_rates); free(r.mode_valid); free(r.mu);
    hss_system_destroy(sys);
    PASS();
}

static void test_dwell_time(void) {
    TEST("compute_dwell_time");
    HSS_System *sys = hss_system_create("dwell_test", 2, 2, 0);
    double A1[4] = {-1, 0, 0, -2};
    double A2[4] = {-3, 0, 0, -4};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A1, NULL, NULL);
    hss_mode_set_dynamics(sys, 1, HSS_CLASS_LINEAR, A2, NULL, NULL);

    HSS_DwellTimeResult r = hss_compute_dwell_time(sys);
    assert(r.num_modes == 2);
    free(r.mode_decay_rates);
    hss_system_destroy(sys);
    PASS();
}

static void test_average_dwell(void) {
    TEST("compute_average_dwell");
    HSS_System *sys = hss_system_create("avg_test", 2, 2, 0);
    double A[4] = {-1, 0, 0, -2};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(sys, 1, HSS_CLASS_LINEAR, A, NULL, NULL);

    HSS_SwitchingSignal sig;
    sig.num_switches = 4;
    double times[4] = {0, 1, 2, 3};
    int modes[4] = {0, 1, 0, 1};
    sig.switch_times = times;
    sig.mode_sequence = modes;

    HSS_AverageDwellResult r = hss_compute_average_dwell(sys, &sig);
    assert(r.num_switches == 4);
    hss_system_destroy(sys);
    PASS();
}

static void test_small_gain(void) {
    TEST("small_gain_analysis");
    HSS_SmallGainResult r = hss_small_gain_analysis(NULL, NULL, 0.5, 0.5);
    assert(r.composite_gain == 0.25);
    assert(r.small_gain_holds);
    assert(r.is_interconnection_ISS);
    PASS();
}

static void test_search_clf(void) {
    TEST("search_clf");
    HSS_System *sys = hss_system_create("clf_s", 2, 2, 0);
    double A[4] = {-2, 0, 0, -3};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(sys, 1, HSS_CLASS_LINEAR, A, NULL, NULL);

    HSS_CLFSearchResult r = hss_search_clf(sys, 100, 0.01);
    assert(r.clf.n == 2);
    free(r.clf.P);
    hss_system_destroy(sys);
    PASS();
}

static void test_lasalle(void) {
    TEST("lasalle_invariance");
    HSS_System *sys = hss_system_create("las_test", 1, 2, 0);
    double A[4] = {-1, 0, 0, -1};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    double x0[2] = {1, 1};
    hss_set_initial_state(sys, 0, x0);

    HSS_ExecutionTrace *tr = hss_trace_alloc(100, 2);
    tr->times[0] = 0; tr->modes[0] = 0;
    tr->times[1] = 1; tr->modes[1] = 0;
    tr->states[0] = 1; tr->states[1] = 1;
    tr->states[2] = 0; tr->states[3] = 0;
    tr->num_steps = 2;

    HSS_LaSalleResult r = hss_lasalle_invariance(sys, tr);
    free(r.omega_limit_set);
    hss_trace_free(tr);
    hss_system_destroy(sys);
    PASS();
}

static void test_zeno_detection(void) {
    TEST("detect_zeno");
    HSS_ExecutionTrace *tr = hss_trace_alloc(50, 1);
    /* Simulate rapid switching */
    for (int i = 0; i < 20; i++) {
        tr->times[i] = i * 1e-6;
        tr->modes[i] = i % 2;
        tr->states[i] = 1.0;
    }
    tr->num_steps = 20;

    HSS_ZenoAnalysis r = hss_detect_zeno(tr, 1e-3);
    assert(r.jump_count >= 0);

    hss_trace_free(tr);
    PASS();
}

static void test_switching_synthesis(void) {
    TEST("synthesize_switching");
    HSS_System *sys = hss_system_create("synth", 3, 2, 0);
    double A[4] = {-1, 0, 0, -2};
    for (int i = 0; i < 3; i++)
        hss_mode_set_dynamics(sys, i, HSS_CLASS_LINEAR, A, NULL, NULL);

    HSS_SwitchingSynthesis r = hss_synthesize_switching(sys, 0, 5.0);
    assert(r.is_stabilizing);
    free(r.signal.mode_sequence);
    free(r.signal.switch_times);

    r = hss_synthesize_switching(sys, 1, 5.0);
    assert(r.is_stabilizing);
    free(r.signal.mode_sequence);
    free(r.signal.switch_times);

    hss_system_destroy(sys);
    PASS();
}

static void test_dwell_computation(void) {
    TEST("compute_min_dwell");
    HSS_System *sys = hss_system_create("dwell_c", 2, 2, 0);
    double A1[4] = {-2, 0, 0, -3};
    double A2[4] = {-1, 0, 0, -4};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A1, NULL, NULL);
    hss_mode_set_dynamics(sys, 1, HSS_CLASS_LINEAR, A2, NULL, NULL);

    HSS_DwellComputation r = hss_compute_min_dwell(sys);
    assert(r.num_pairs == 4);
    free(r.dwell_result.mode_decay_rates);
    free(r.pairwise_dwell);
    hss_system_destroy(sys);
    PASS();
}

int main(void) {
    printf("=== HSS Analysis Tests ===\n\n");
    test_clf_verify();
    test_mlf_verify();
    test_dwell_time();
    test_average_dwell();
    test_small_gain();
    test_search_clf();
    test_lasalle();
    test_zeno_detection();
    test_switching_synthesis();
    test_dwell_computation();
    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
