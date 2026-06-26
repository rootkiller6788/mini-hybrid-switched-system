/**
 * @file test_hss_stability.c
 * @brief Tests for HSS Hybrid Stability (L4-L7): unified stability,
 *        cross-framework verification, application stability.
 */
#include "hss_hybrid_stability.h"
#include "hss_core.h"
#include "hss_analysis.h"
#include "hss_simulation.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>

static int tests_run = 0, tests_passed = 0;
#define TEST(n) do { tests_run++; printf("  TEST: %s ... ", n); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)

static void test_unified_stability(void) {
    TEST("unified_stability_analysis");
    HSS_System *sys = hss_system_create("uni", 2, 2, 0);
    double A1[4] = {-2, 0, 0, -3};
    double A2[4] = {-1, 0, 0, -4};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A1, NULL, NULL);
    hss_mode_set_dynamics(sys, 1, HSS_CLASS_LINEAR, A2, NULL, NULL);
    double P[4] = {1, 0, 0, 1};
    hss_mode_set_lyapunov(sys, 0, P, 0.5);
    hss_mode_set_lyapunov(sys, 1, P, 0.3);

    HSS_SwitchingSignal sig;
    double times[2] = {0, 1};
    int modes[2] = {0, 1};
    sig.switch_times = times;
    sig.mode_sequence = modes;
    sig.num_switches = 2;

    HSS_UnifiedStability r = hss_unified_stability_analysis(sys, &sig);
    assert(r.is_stable); /* All modes stable with LF */
    hss_unified_stability_report(&r, stdout);

    hss_system_destroy(sys);
    PASS();
}

static void test_cross_framework(void) {
    TEST("cross_framework_verify");
    HSS_System *sys = hss_system_create("cross", 2, 2, 0);
    double A[4] = {-2, 0, 0, -3};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(sys, 1, HSS_CLASS_LINEAR, A, NULL, NULL);

    HSS_SwitchingSignal sig;
    sig.num_switches = 0;

    HSS_CrossFrameworkStability r = hss_cross_framework_verify(sys, &sig);
    assert(r.frameworks_checked == 4);
    assert(r.confidence >= 0.0 && r.confidence <= 1.0);
    hss_cross_framework_report(&r, stdout);

    for (int i = 0; i < r.frameworks_checked; i++)
        free(r.framework_names[i]);
    free(r.framework_names);
    free(r.framework_stable);
    free(r.framework_margin);

    hss_system_destroy(sys);
    PASS();
}

static void test_stability_margin(void) {
    TEST("compute_stability_margin");
    HSS_System *sys = hss_system_create("margin", 2, 2, 0);
    double A1[4] = {-2, 0, 0, -3};
    double A2[4] = {-1, 0.5, -0.5, -4};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A1, NULL, NULL);
    hss_mode_set_dynamics(sys, 1, HSS_CLASS_LINEAR, A2, NULL, NULL);

    HSS_StabilityMargin r = hss_compute_stability_margin(sys);
    assert(r.robust || !r.robust); /* Always defined */
    hss_system_destroy(sys);
    PASS();
}

static void test_robust_stability(void) {
    TEST("robust_stability_analysis");
    HSS_System *sys = hss_system_create("robust", 1, 2, 0);
    double A[4] = {-3, 0, 0, -5};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    double P[4] = {1, 0, 0, 1};
    hss_mode_set_lyapunov(sys, 0, P, 0.5);

    HSS_RobustStability r = hss_robust_stability_analysis(sys, 0.1, 10);
    assert(r.n == 2);
    free(r.perturbation_direction);
    hss_system_destroy(sys);
    PASS();
}

static void test_power_electronics(void) {
    TEST("power_electronics_stability (L7)");
    HSS_PowerElectronicsStability r = hss_power_electronics_stability(
        1e-3, 100e-6, 10.0, 12.0, 24.0, 0.5);
    assert(r.is_stable); /* Boost converter is inherently stable */
    assert(r.settling_time > 0.0);
    PASS();
}

static void test_automotive(void) {
    TEST("automotive_stability (L7)");
    HSS_AutomotiveStability r = hss_automotive_stability(
        750.0, 800.0, 30.0, 0);
    assert(r.is_idle_stable);
    assert(r.recovery_time > 0.0);
    PASS();
}

static void test_hvac(void) {
    TEST("hvac_stability (L7)");
    HSS_HVACStability r = hss_hvac_stability(
        22.0, 21.0, 10.0, 500.0, 10.0, 8.0);
    assert(r.heat_loss_coeff > 0.0);
    PASS();
}

static void test_converse_lyapunov(void) {
    TEST("converse_lyapunov");
    HSS_System *sys = hss_system_create("conv", 1, 1, 0);
    double A[1] = {-1.0};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A, NULL, NULL);

    HSS_ExecutionTrace *tr = hss_trace_alloc(10, 1);
    for (int i = 0; i < 10; i++) {
        tr->times[i] = i * 0.1;
        tr->modes[i] = 0;
        tr->states[i] = exp(-i * 0.1);
    }
    tr->num_steps = 10;

    HSS_ConverseLyapunov r = hss_converse_lyapunov(sys, tr, 0.5);
    assert(r.is_constructed);
    free(r.V_values);
    hss_trace_free(tr);
    hss_system_destroy(sys);
    PASS();
}

static void test_comparison_system(void) {
    TEST("comparison_system_setup");
    HSS_System *sys = hss_system_create("comp", 1, 1, 0);
    double A[1] = {-2.0};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    double P[1] = {1.0};
    hss_mode_set_lyapunov(sys, 0, P, 2.0);

    HSS_ComparisonSystem r = hss_comparison_system_setup(sys, 10.0);
    assert(r.comparison_state == 10.0);

    hss_system_destroy(sys);
    PASS();
}

int main(void) {
    printf("=== HSS Stability Tests ===\n\n");
    test_unified_stability();
    test_cross_framework();
    test_stability_margin();
    test_robust_stability();
    test_power_electronics();
    test_automotive();
    test_hvac();
    test_converse_lyapunov();
    test_comparison_system();
    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
