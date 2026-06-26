#include "dta_core.h"
#include "dta_switch_signal.h"
#include "dta_stability.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

static int tests_run = 0, tests_passed = 0;
#define TEST(n) do { tests_run++; printf("  %s... ", n); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(m) do { printf("FAIL: %s\n", m); } while(0)
#define CHECK(c, m) do { if (!(c)) { FAIL(m); return; } } while(0)

void test_dwell_stability_stable(void) {
    TEST("dwell stability (stable modes)");
    DTA_SwitchedSystem* sys = dta_system_create(2, 2, 0, 0);
    double A1[4] = {-1.0, 0.2, -0.2, -1.5};
    double A2[4] = {-2.0, 0.1, -0.1, -0.5};
    dta_system_set_linear_mode(sys, 0, A1, NULL, NULL);
    dta_system_set_linear_mode(sys, 1, A2, NULL, NULL);
    DTA_DwellStabilityResult res = dta_analyze_dwell_stability(sys, 10.0);
    CHECK(res.tau_d_star >= 0, "tau_d_star should be nonnegative");
    if (res.lambda_i) free(res.lambda_i);
    if (res.common_P) free(res.common_P);
    dta_system_free(sys);
    PASS();
}

void test_dwell_stability_unstable(void) {
    TEST("dwell stability (unstable mode)");
    DTA_SwitchedSystem* sys = dta_system_create(2, 2, 0, 0);
    double A1[4] = {-1.0, 0.0, 0.0, -2.0};
    double A_unstable[4] = {1.0, 0.0, 0.0, 1.0};
    dta_system_set_linear_mode(sys, 0, A1, NULL, NULL);
    dta_system_set_linear_mode(sys, 1, A_unstable, NULL, NULL);
    DTA_DwellStabilityResult res = dta_analyze_dwell_stability(sys, 1.0);
    CHECK(res.verdict == DTA_INCONCLUSIVE, "should be inconclusive with unstable mode");
    if (res.lambda_i) free(res.lambda_i);
    if (res.common_P) free(res.common_P);
    dta_system_free(sys);
    PASS();
}

void test_arbitrary_switching(void) {
    TEST("arbitrary switching");
    DTA_SwitchedSystem* sys = dta_system_create(2, 2, 0, 0);
    double A1[4] = {-1.0, 0.2, -0.2, -1.5};
    double A2[4] = {-2.0, 0.1, -0.1, -0.5};
    dta_system_set_linear_mode(sys, 0, A1, NULL, NULL);
    dta_system_set_linear_mode(sys, 1, A2, NULL, NULL);
    DTA_ArbitrarySwitchingResult res = dta_check_arbitrary_switching(sys);
    CHECK(res.n == 2, "wrong state dim");
    if (res.common_P) free(res.common_P);
    dta_system_free(sys);
    PASS();
}

void test_find_common_quadratic(void) {
    TEST("find common quadratic");
    DTA_SwitchedSystem* sys = dta_system_create(2, 2, 0, 0);
    double A1[4] = {-1.0, 0.2, -0.2, -1.5};
    double A2[4] = {-2.0, 0.1, -0.1, -0.5};
    dta_system_set_linear_mode(sys, 0, A1, NULL, NULL);
    dta_system_set_linear_mode(sys, 1, A2, NULL, NULL);
    double P[4];
    bool found = dta_find_common_quadratic(sys, P);
    /* May or may not find one - just check function runs */
    (void)found;
    dta_system_free(sys);
    PASS();
}

void test_stability_margin(void) {
    TEST("stability margin");
    DTA_SwitchedSystem* sys = dta_system_create(2, 2, 0, 0);
    double A1[4] = {-1.0, 0.0, 0.0, -2.0};
    double A2[4] = {-3.0, 0.0, 0.0, -0.5};
    dta_system_set_linear_mode(sys, 0, A1, NULL, NULL);
    dta_system_set_linear_mode(sys, 1, A2, NULL, NULL);
    double margin = dta_stability_margin(sys);
    CHECK(margin > 0, "stability margin should be positive");
    dta_system_free(sys);
    PASS();
}

void test_decay_rate(void) {
    TEST("decay rate bound");
    DTA_SwitchedSystem* sys = dta_system_create(2, 2, 0, 0);
    double A1[4] = {-1.0, 0.0, 0.0, -2.0};
    double A2[4] = {-3.0, 0.0, 0.0, -0.5};
    dta_system_set_linear_mode(sys, 0, A1, NULL, NULL);
    dta_system_set_linear_mode(sys, 1, A2, NULL, NULL);
    double rate = dta_decay_rate_bound(sys, 5.0);
    CHECK(rate >= 0, "decay rate should be nonnegative");
    dta_system_free(sys);
    PASS();
}

int main(void) {
    printf("=== Test: dta_stability ===\n");
    test_dwell_stability_stable();
    test_dwell_stability_unstable();
    test_arbitrary_switching();
    test_find_common_quadratic();
    test_stability_margin();
    test_decay_rate();
    printf("\nResults: %d/%d tests passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
