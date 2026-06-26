#include "dta_core.h"
#include "dta_switch_signal.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

/* test_signal.c - Tests for switching signal generation and analysis */

static int tests_run = 0, tests_passed = 0;
#define TEST(n) do { tests_run++; printf("  %s... ", n); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(m) do { printf("FAIL: %s\n", m); } while(0)
#define CHECK(c, m) do { if (!(c)) { FAIL(m); return; } } while(0)

void test_periodic_signal(void) {
    TEST("periodic signal");
    int pattern[] = {0, 1};
    DTA_SwitchingSignal* sig = dta_signal_periodic(0.0, 5.0, pattern, 2, 1.0);
    CHECK(sig != NULL, "null signal");
    CHECK(sig->n_switches >= 4, "too few switches");
    CHECK(sig->mode_sequence[0] == 0, "wrong first mode");
    CHECK(sig->mode_sequence[1] == 1, "wrong second mode");
    dta_signal_free(sig);
    PASS();
}

void test_constant_dwell(void) {
    TEST("constant dwell signal");
    DTA_SwitchingSignal* sig = dta_signal_constant_dwell(0.0, 10.0, 2, 2.0);
    CHECK(sig != NULL, "null signal");
    double min_d = dta_signal_min_dwell(sig);
    CHECK(min_d >= 2.0 - 1e-6, "dwell not respected");
    dta_signal_free(sig);
    PASS();
}

void test_random_dwell(void) {
    TEST("random dwell signal");
    srand(42);
    DTA_SwitchingSignal* sig = dta_signal_random_dwell(0.0, 10.0, 3, 1.0);
    CHECK(sig != NULL, "null signal");
    CHECK(sig->n_switches >= 1, "no switches");
    dta_signal_free(sig);
    PASS();
}

void test_signal_statistics(void) {
    TEST("signal statistics");
    int pat[] = {0, 1, 0, 1};
    DTA_SwitchingSignal* sig = dta_signal_periodic(0.0, 8.0, pat, 2, 1.0);
    DTA_SignalStatistics stats = dta_signal_statistics(sig, 2);
    CHECK(stats.total_switches > 0, "zero switches");
    CHECK(stats.duration > 0, "zero duration");
    CHECK(stats.switch_frequency > 0, "zero frequency");
    dta_signal_free(sig);
    PASS();
}

void test_active_mode(void) {
    TEST("active mode query");
    int pat[] = {0, 1};
    DTA_SwitchingSignal* sig = dta_signal_periodic(0.0, 5.0, pat, 2, 1.0);
    CHECK(dta_signal_active_mode(sig, 0.5) == 0, "wrong mode at 0.5");
    CHECK(dta_signal_active_mode(sig, 1.5) == 1, "wrong mode at 1.5");
    dta_signal_free(sig);
    PASS();
}

void test_signal_validate_dwell(void) {
    TEST("validate dwell");
    DTA_SwitchingSignal* sig = dta_signal_constant_dwell(0.0, 10.0, 2, 2.0);
    CHECK(dta_signal_validate_dwell(sig, 1.9), "should pass at 1.9");
    CHECK(dta_signal_validate_dwell(sig, 2.0), "should pass at 2.0");
    CHECK(!dta_signal_validate_dwell(sig, 2.1), "should fail at 2.1");
    dta_signal_free(sig);
    PASS();
}

void test_adt_bound(void) {
    TEST("ADT bound check");
    DTA_SwitchingSignal* sig = dta_signal_constant_dwell(0.0, 10.0, 2, 2.0);
    bool ok = dta_signal_check_adt_bound(sig, 2.0, 1.0, 10.0, 0.0);
    CHECK(ok, "ADT bound should hold");
    dta_signal_free(sig);
    PASS();
}

void test_signal_slice(void) {
    TEST("signal slice");
    int pat[] = {0, 1, 0, 1};
    DTA_SwitchingSignal* sig = dta_signal_periodic(0.0, 10.0, pat, 2, 1.0);
    DTA_SwitchingSignal* sub = dta_signal_slice(sig, 2.0, 5.0);
    CHECK(sub != NULL, "null slice");
    CHECK(sub->n_switches > 0, "empty slice");
    dta_signal_free(sig);
    dta_signal_free(sub);
    PASS();
}

int main(void) {
    printf("=== Test: dta_switch_signal ===\n");
    test_periodic_signal();
    test_constant_dwell();
    test_random_dwell();
    test_signal_statistics();
    test_active_mode();
    test_signal_validate_dwell();
    test_adt_bound();
    test_signal_slice();
    printf("\nResults: %d/%d tests passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
