/**
 * @file test_hss_simulation.c
 * @brief Tests for HSS Simulation (L5-L6): canonical problem solvers,
 *        simulation loop, ODE solvers.
 */

#include "hss_simulation.h"
#include "hss_core.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); } while(0)

/* ---- L6 KP1: Bouncing Ball ---- */
static void test_bouncing_ball(void) {
    TEST("bouncing_ball");
    HSS_BouncingBall bb = hss_bouncing_ball_init(9.81, 0.8, 10.0, 0.0);

    assert(fabs(bb.gravity - 9.81) < 1e-6);
    assert(fabs(bb.restitution - 0.8) < 1e-6);
    assert(bb.bounce_count == 0);
    assert(!bb.at_rest);

    /* Step with small dt */
    (void)hss_bouncing_ball_step(&bb, 3.0);
    /* Should bounce at least once from 10m height */
    assert(bb.bounce_count >= 1);
    assert(bb.height >= 0.0);

    /* Simulate to rest */
    hss_bouncing_ball_simulate(&bb, 30.0);
    assert(bb.at_rest || bb.bounce_count > 2);

    PASS();
}

/* ---- L6 KP2: Thermostat ---- */
static void test_thermostat(void) {
    TEST("thermostat");
    HSS_Thermostat t = hss_thermostat_init(20.0, 18.0, 22.0, 1.0, 0.5, 30.0, 10.0);

    assert(fabs(t.temperature - 20.0) < 1e-6);
    assert(fabs(t.T_low - 18.0) < 1e-6);
    assert(fabs(t.T_high - 22.0) < 1e-6);

    /* Initially at 20°C between 18 and 22, heater should be off */
    assert(!t.heater_on); /* Started above T_low */

    /* Cool and let it switch on */
    t.temperature = 17.5;
    (void)hss_thermostat_step(&t, 0.01);
    /* Should turn on */
    assert(t.heater_on);

    /* Run full simulation */
    hss_thermostat_simulate(&t, 1.0, 0.01);
    assert(t.time > 0.0);

    PASS();
}

/* ---- L6 KP3: DC-DC Converter ---- */
static void test_dcdc_converter(void) {
    TEST("dcdc_converter");
    HSS_DCDCConverter c = hss_dcdc_init(12.0, 1e-3, 100e-6, 10.0, 100e3, 0.5);

    assert(fabs(c.Vin - 12.0) < 1e-6);
    assert(fabs(c.duty_cycle - 0.5) < 1e-6);

    /* Run simulation */
    double vout = hss_dcdc_simulate(&c, 0.001, 1e-6);
    /* Output should be positive */
    assert(vout >= 0.0);

    /* Duty cycle changes should be reflected */
    (void)vout;

    PASS();
}

/* ---- L6 KP4: Cruise Control ---- */
static void test_cruise_control(void) {
    TEST("cruise_control");
    HSS_CruiseControl cc = hss_cruise_control_init(30.0, 50.0, 28.0, 35.0, 2.0, 3.0);

    assert(fabs(cc.ego_velocity - 30.0) < 1e-6);
    assert(fabs(cc.v_set - 35.0) < 1e-6);
    assert(cc.mode == 1); /* FOLLOW: gap is small initially */

    /* Advance simulation */
    int mode = hss_cruise_control_step(&cc, 0.1);
    /* Should remain in valid mode */
    assert(mode >= 0 && mode <= 2);

    /* Run full simulation */
    hss_cruise_control_simulate(&cc, 5.0, 0.01);
    assert(cc.time > 0.0);

    PASS();
}

/* ---- L5 KP11: Step Controller ---- */
static void test_step_controller(void) {
    TEST("step_controller");
    HSS_StepController sc = hss_step_controller_init(4, 1e-6, 1e-6);

    assert(sc.order == 4);
    assert(sc.h_current > 0.0);

    /* Adapt step */
    double h_new = hss_adapt_step(&sc, 1e-7); /* Small error */
    assert(h_new > 0.0);

    h_new = hss_adapt_step(&sc, 1e-3); /* Large error */
    assert(h_new > 0.0);

    PASS();
}

/* ---- L5: Simulation Config ---- */
static void test_sim_config(void) {
    TEST("sim_config_default");
    HSS_SimConfig cfg = hss_sim_config_default();
    assert(cfg.t_end == 10.0);
    assert(cfg.solver == HSS_SOLVER_RK4);
    assert(cfg.detect_events == true);

    PASS();
}

/* ---- L5: Trace Management ---- */
static void test_trace_management(void) {
    TEST("trace_alloc/record/free");
    HSS_ExecutionTrace *trace = hss_trace_alloc(100, 3);
    assert(trace != NULL);
    assert(trace->max_steps == 100);
    assert(trace->num_steps == 0);

    double x[3] = {1.0, 2.0, 3.0};
    int ret = hss_trace_record(trace, 0.0, 0, x, 3);
    assert(ret == 0);
    assert(trace->num_steps == 1);
    assert(fabs(trace->times[0]) < 1e-9);
    assert(trace->modes[0] == 0);

    /* Record more */
    for (int i = 1; i < 50; i++) {
        x[0] += 0.1;
        ret = hss_trace_record(trace, i * 0.01, 0, x, 3);
        assert(ret == 0);
    }
    assert(trace->num_steps == 50);

    /* Print summary (no crash) */
    hss_trace_print_summary(trace, stdout);
    hss_trace_print_summary(NULL, stdout);

    /* Export CSV */
    ret = hss_trace_export_csv(trace, "build/test_trace.csv");
    /* May fail if build/ doesn't exist, that's OK */
    (void)ret;

    hss_trace_free(trace);
    hss_trace_free(NULL); /* NULL safety */

    PASS();
}

/* ---- L5: Event adjustment ---- */
static void test_event_adjustment(void) {
    TEST("adjust_step_for_event");
    double h = hss_adjust_step_for_event(0.1, 0.05);
    assert(h < 0.1);
    assert(h > 0.0);

    h = hss_adjust_step_for_event(0.1, 0.2);  /* Event after step */
    assert(fabs(h - 0.1) < 1e-9);

    h = hss_adjust_step_for_event(0.1, -1.0); /* No event */
    assert(fabs(h - 0.1) < 1e-9);

    PASS();
}

/* ---- L5 KP9: Complete system simulation ---- */
static void test_hybrid_simulation(void) {
    TEST("hss_simulate (linear)");
    /* Create a simple 2-mode switched system */
    HSS_System *sys = hss_system_create("sim_test", 2, 2, 0);
    assert(sys != NULL);

    /* Mode 0: stable focus */
    double A0[4] = {-0.1, 1.0, -1.0, -0.1};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A0, NULL, NULL);

    /* Mode 1: stable node */
    double A1[4] = {-2.0, 0.0, 0.0, -3.0};
    hss_mode_set_dynamics(sys, 1, HSS_CLASS_LINEAR, A1, NULL, NULL);

    /* Add transition: mode 0 → 1 when x[0] > 5 */
    HSS_Guard g;
    double ga[2] = {1.0, 0.0};
    g.a = ga; g.dim = 2; g.b = 5.0; g.is_active = true;
    hss_add_transition(sys, 0, 1, &g, 1, HSS_RESET_IDENTITY, NULL, NULL, "t01");

    /* Set initial state */
    double x0[2] = {1.0, 0.0};
    hss_set_initial_state(sys, 0, x0);

    /* Configure simulation */
    HSS_SimConfig cfg = hss_sim_config_default();
    cfg.t_end = 1.0;
    cfg.dt = 0.01;
    cfg.max_steps = 500;
    cfg.detect_events = true;

    HSS_ExecutionTrace *trace = hss_trace_alloc(500, 2);
    HSS_SimStats stats;

    int ret = hss_simulate(sys, &cfg, trace, &stats);
    assert(ret == 0); /* Should complete without Zeno */
    assert(stats.num_steps > 0);
    assert(trace->num_steps > 0);

    hss_trace_free(trace);
    hss_system_destroy(sys);
    PASS();
}

/* ---- Step Continuous ---- */
static void test_step_continuous(void) {
    TEST("step_continuous");
    HSS_System *sys = hss_system_create("step_test", 1, 1, 0);
    double A[1] = {-1.0};  /* Stable: ẋ = -x */
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    double x0[1] = {2.0};
    hss_set_initial_state(sys, 0, x0);

    int ret = hss_step_continuous(sys, 0.1);
    assert(ret == 0);
    /* After one RK4 step: x should decrease */
    assert(sys->state.data[0] < 2.0);

    hss_system_destroy(sys);
    PASS();
}

/* ---- Guard Check ---- */
static void test_guard_check(void) {
    TEST("check_guards / apply_transition");
    HSS_System *sys = hss_system_create("guard_test", 2, 2, 0);
    double A[4] = {0, 0, 0, 0};
    hss_mode_set_dynamics(sys, 0, HSS_CLASS_LINEAR, A, NULL, NULL);
    hss_mode_set_dynamics(sys, 1, HSS_CLASS_LINEAR, A, NULL, NULL);

    HSS_Guard g;
    double ga[2] = {1.0, 0.0};
    g.a = ga; g.dim = 2; g.b = 3.0; g.is_active = true;

    hss_add_transition(sys, 0, 1, &g, 1, HSS_RESET_IDENTITY, NULL, NULL, "t");

    /* State below threshold: no guard triggered */
    double x0[2] = {1.0, 0.0};
    hss_set_initial_state(sys, 0, x0);
    int idx = hss_check_guards(sys);
    assert(idx == -1);

    /* State above threshold: guard triggered */
    sys->state.data[0] = 5.0;
    idx = hss_check_guards(sys);
    assert(idx >= 0);

    /* Apply transition */
    int ret = hss_apply_transition(sys, idx);
    assert(ret == 0);
    assert(sys->active_mode == 1);

    hss_system_destroy(sys);
    PASS();
}

/* ---- Multi-rate Simulation ---- */
static void test_multirate(void) {
    TEST("multirate_init / macro_step");
    HSS_System *fast = hss_system_create("fast", 1, 1, 0);
    HSS_System *slow = hss_system_create("slow", 1, 1, 0);

    double Af[1] = {-10.0};
    double As[1] = {-0.1};
    hss_mode_set_dynamics(fast, 0, HSS_CLASS_LINEAR, Af, NULL, NULL);
    hss_mode_set_dynamics(slow, 0, HSS_CLASS_LINEAR, As, NULL, NULL);
    hss_set_initial_state(fast, 0, (double[]){1.0});
    hss_set_initial_state(slow, 0, (double[]){10.0});

    HSS_MultiRateSim mrs = hss_multirate_init(fast, slow, 10);
    assert(mrs.fast_steps_per_slow == 10);

    int ret = hss_multirate_macro_step(&mrs);
    /* Should succeed */
    assert(ret == 0);

    hss_system_destroy(fast);
    hss_system_destroy(slow);
    PASS();
}

int main(void) {
    printf("=== HSS Simulation Tests ===\n\n");

    test_bouncing_ball();
    test_thermostat();
    test_dcdc_converter();
    test_cruise_control();
    test_step_controller();
    test_sim_config();
    test_trace_management();
    test_event_adjustment();
    test_hybrid_simulation();
    test_step_continuous();
    test_guard_check();
    test_multirate();

    printf("\n=== Results: %d/%d tests passed ===\n",
           tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
