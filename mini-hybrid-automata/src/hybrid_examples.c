/**
 * @file hybrid_examples.c
 * @brief Canonical hybrid automaton examples (KP49-KP55, L6-L7)
 *
 * Standard benchmark examples used in hybrid systems research:
 * bouncing ball, thermostat, two-tank system, train-gate controller,
 * DC-DC converter, engine air-fuel control, robot obstacle avoidance.
 *
 * Each example constructs a complete hybrid automaton with modes,
 * transitions, guards, resets, invariants, and flows.
 *
 * Reference:
 *   Johansson, "Piecewise Linear Control Systems" (2003)
 *   Alur et al., "Algorithmic Analysis of Hybrid Systems" (1995)
 *   Lygeros et al., "Lecture Notes on Hybrid Systems" (2004)
 */

#include "hybrid_automaton.h"
#include "hybrid_execution.h"
#include "hybrid_examples.h"
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * KP49: Bouncing Ball
 * ========================================================================== */

/**
 * @brief KP49: Create the classical bouncing ball hybrid automaton.
 *
 * A ball dropped from height h₀ bounces on the ground (height = 0).
 * Each bounce dissipates energy — velocity reverses with restitution
 * coefficient c ∈ (0, 1).
 *
 * State: x = (h, v) — height above ground, vertical velocity
 *
 * Modes: single mode "Falling"
 *   Flow:    ḣ = v,  v̇ = -g   (gravity)
 *   Inv:     h ≥ 0             (above ground)
 *   Transition: when h = 0 and v < 0 (hitting ground)
 *     Reset: h' = h, v' = -c·v  (bounce)
 *
 * This system exhibits Zeno behavior: infinitely many bounces in
 * finite time. Total time: T_∞ = v₀/g + 2v₀/(g(1-c))
 *
 * @param h0 Initial height (must be > 0)
 * @param v0 Initial velocity (> 0 for upward)
 * @param c  Restitution coefficient (0 < c < 1)
 * @param g  Gravitational acceleration (default 9.81)
 * @return   Hybrid automaton of the bouncing ball
 *
 * Complexity: O(1) construction
 */
HybridAutomaton* example_bouncing_ball(double h0, double v0, double c, double g)
{
    HybridAutomaton *ha = hybrid_automaton_create("BouncingBall", 2);
    if (!ha) return NULL;

    /* Set variables: x[0] = height h, x[1] = velocity v */
    hybrid_automaton_set_variable(ha, 0, "h", HAVAR_REAL_BOUNDED, 0.0, INFINITY);
    hybrid_automaton_set_variable(ha, 1, "v", HAVAR_REAL, -INFINITY, INFINITY);

    /* Mode 0: Falling */
    int falling = hybrid_automaton_add_mode(ha, "Falling", HAMODE_INITIAL);
    if (falling < 0) { hybrid_automaton_destroy(ha); return NULL; }

    /* Flow: ḣ = v, v̇ = -g
       A = [[0, 1], [0, 0]], b = [0, -g] */
    double A[4] = {0, 1, 0, 0};
    double b[2] = {0, -g};
    hybrid_flow_set(ha, falling, HAFLOW_AFFINE, A, b);

    /* Invariant: h ≥ 0 → -h ≤ 0 → H=[-1,0], k=0 */
    double H[2] = {-1, 0};
    double k = 0;
    hybrid_invariant_set(ha, falling, 1, H, &k);

    /* Transition: bounce (self-loop)
       Guard: h ≤ ε (near ground) and v ≤ 0 (falling)
       Use constraints: [1,0]*x ≤ ε (h ≤ ε) and [0,1]*x ≤ 0 (v ≤ 0)
       With ε=0.01 for numerical tolerance in discrete-time simulation. */
    int bounce = hybrid_automaton_add_transition(ha, falling, falling,
                                                   HATRIG_AUTONOMOUS, "bounce");
    if (bounce < 0) { hybrid_automaton_destroy(ha); return NULL; }

    /* ε = 0.1: guard triggers when h ≤ 0.1m, accommodating discrete-time steps */
    double eps = 0.1;
    double guard_A[4] = {1, 0,   0, 1};
    double guard_b[2] = {eps, 0};
    hybrid_guard_set(ha, bounce, 2, guard_A, guard_b);

    /* Reset: h' = h, v' = -c·v
       R = [[1, 0], [0, -c]], r = [0, 0] */
    double R[4] = {1, 0, 0, -c};
    double r[2] = {0, 0};
    hybrid_reset_set(ha, bounce, HARESET_AFFINE, R, r);

    /* Initial condition: h = h₀, v = v₀ */
    double x0[2] = {h0, v0};
    hybrid_init_set(ha, falling, x0);

    return ha;
}

/* ==========================================================================
 * KP50: Thermostat
 * ========================================================================== */

/**
 * @brief KP50: Create the thermostat hybrid automaton.
 *
 * A thermostat controls heating to maintain temperature within
 * a target range [T_lo, T_hi].
 *
 * State: x = (T) — room temperature
 *
 * Modes:
 *   ON:  Heater active,   Ṫ = K_h·(T_env - T) + K_on  (heating)
 *   OFF: Heater inactive, Ṫ = K_h·(T_env - T)           (cooling)
 *
 * Transitions:
 *   ON → OFF: when T ≥ T_hi  (too warm)
 *   OFF → ON: when T ≤ T_lo  (too cold)
 *   Guard: T boundary crossing
 *   Reset: identity (T continuous)
 *
 * @param T_lo   Low temperature threshold
 * @param T_hi   High temperature threshold
 * @param T_env  Environmental temperature
 * @param K_h    Heat transfer coefficient
 * @param K_on   Heating rate
 * @param T_init Initial temperature
 * @param start_on Initial mode is ON if true
 * @return       Hybrid automaton of the thermostat
 */
HybridAutomaton* example_thermostat(double T_lo, double T_hi, double T_env,
                                     double K_h, double K_on,
                                     double T_init, bool start_on)
{
    HybridAutomaton *ha = hybrid_automaton_create("Thermostat", 1);
    if (!ha) return NULL;

    hybrid_automaton_set_variable(ha, 0, "T", HAVAR_REAL, -INFINITY, INFINITY);

    /* Mode 0: ON (heating) */
    int mode_on = hybrid_automaton_add_mode(ha, "ON",
                                              start_on ? HAMODE_INITIAL : HAMODE_NORMAL);
    /* Affine form: Ṫ = -K_h·T + (K_h·T_env + K_on) */
    double A[1] = {-K_h};
    double b_on_vec[1] = {K_h * T_env + K_on};
    hybrid_flow_set(ha, mode_on, HAFLOW_AFFINE, A, b_on_vec);

    /* Invariant: T ≤ T_hi (so ON mode is bounded above) */
    double H_hi[1] = {1};
    double k_hi = T_hi;
    hybrid_invariant_set(ha, mode_on, 1, H_hi, &k_hi);

    /* Mode 1: OFF (cooling) */
    int mode_off = hybrid_automaton_add_mode(ha, "OFF",
                                               start_on ? HAMODE_NORMAL : HAMODE_INITIAL);
    double b_off_vec[1] = {K_h * T_env};
    hybrid_flow_set(ha, mode_off, HAFLOW_AFFINE, A, b_off_vec);

    /* Invariant: T ≥ T_lo */
    double H_lo[1] = {-1};
    double k_lo = -T_lo;
    hybrid_invariant_set(ha, mode_off, 1, H_lo, &k_lo);

    /* Transition: ON → OFF when T ≥ T_hi */
    int on_to_off = hybrid_automaton_add_transition(ha, mode_on, mode_off,
                                                      HATRIG_AUTONOMOUS, "too_hot");
    double guard_hi[1] = {1};
    double guard_b_hi = T_hi;
    hybrid_guard_set(ha, on_to_off, 1, guard_hi, &guard_b_hi);
    hybrid_reset_set(ha, on_to_off, HARESET_IDENTITY, NULL, NULL);

    /* Transition: OFF → ON when T ≤ T_lo */
    int off_to_on = hybrid_automaton_add_transition(ha, mode_off, mode_on,
                                                      HATRIG_AUTONOMOUS, "too_cold");
    double guard_lo[1] = {-1};
    double guard_b_lo = -T_lo;
    hybrid_guard_set(ha, off_to_on, 1, guard_lo, &guard_b_lo);
    hybrid_reset_set(ha, off_to_on, HARESET_IDENTITY, NULL, NULL);

    /* Initial condition */
    double x0[1] = {T_init};
    hybrid_init_set(ha, start_on ? mode_on : mode_off, x0);

    return ha;
}

/* ==========================================================================
 * KP51: Two-Tank System
 * ========================================================================== */

/**
 * @brief KP51: Create the two-tank hybrid automaton.
 *
 * Two connected water tanks with a valve that can be OPEN or CLOSED.
 * When open, water flows from tank 1 to tank 2 at rate proportional
 * to pressure difference.
 *
 * State: x = (h₁, h₂) — water heights
 *
 * Modes:
 *   CLOSED: ḣ₁ = 0, ḣ₂ = -K_out·h₂  (tank 2 drains)
 *   OPEN:   ḣ₁ = -K₁·(h₁-h₂), ḣ₂ = K₁·(h₁-h₂) - K_out·h₂
 *
 * Transitions: controlled (open/close valve)
 *
 * @param K_inter  Inter-tank flow coefficient
 * @param K_out    Output flow coefficient
 * @param h1_init  Initial height tank 1
 * @param h2_init  Initial height tank 2
 * @return         Hybrid automaton
 */
HybridAutomaton* example_two_tank(double K_inter, double K_out,
                                   double h1_init, double h2_init)
{
    HybridAutomaton *ha = hybrid_automaton_create("TwoTank", 2);
    if (!ha) return NULL;

    hybrid_automaton_set_variable(ha, 0, "h1", HAVAR_REAL_BOUNDED, 0.0, INFINITY);
    hybrid_automaton_set_variable(ha, 1, "h2", HAVAR_REAL_BOUNDED, 0.0, INFINITY);

    /* Mode 0: CLOSED */
    int closed = hybrid_automaton_add_mode(ha, "Closed", HAMODE_INITIAL);
    double A_closed[4] = {0, 0, 0, -K_out};
    double b_closed[2] = {0, 0};
    hybrid_flow_set(ha, closed, HAFLOW_LINEAR, A_closed, b_closed);

    /* Invariant: h₁ ≥ 0, h₂ ≥ 0 */
    double H_inv[4] = {-1, 0, 0, -1};
    double k_inv[2] = {0, 0};
    hybrid_invariant_set(ha, closed, 2, H_inv, k_inv);

    /* Mode 1: OPEN */
    int open = hybrid_automaton_add_mode(ha, "Open", HAMODE_NORMAL);
    double A_open[4] = {-K_inter, K_inter, K_inter, -(K_inter + K_out)};
    double b_open[2] = {0, 0};
    hybrid_flow_set(ha, open, HAFLOW_LINEAR, A_open, b_open);
    hybrid_invariant_set(ha, open, 2, H_inv, k_inv);

    /* Transition: CLOSED → OPEN (control) */
    int c2o = hybrid_automaton_add_transition(ha, closed, open,
                                                HATRIG_CONTROLLED, "open_valve");
    hybrid_reset_set(ha, c2o, HARESET_IDENTITY, NULL, NULL);

    /* Transition: OPEN → CLOSED (control) */
    int o2c = hybrid_automaton_add_transition(ha, open, closed,
                                                HATRIG_CONTROLLED, "close_valve");
    hybrid_reset_set(ha, o2c, HARESET_IDENTITY, NULL, NULL);

    double x0[2] = {h1_init, h2_init};
    hybrid_init_set(ha, closed, x0);

    return ha;
}

/* ==========================================================================
 * KP52: Train-Gate Controller (benchmark)
 * ========================================================================== */

/**
 * @brief KP52: Create the train-gate controller hybrid automaton.
 *
 * Classic hybrid verification benchmark: a train approaches a crossing,
 * a gate must be lowered before the train arrives, and raised after it
 * passes. Safety: gate is down whenever train is in crossing.
 *
 * State: x = (train_pos, gate_angle)
 *
 * Modes:
 *   FAR:    train approaching, gate up
 *   NEAR:   train near, gate lowering
 *   IN:     train in crossing, gate down
 *   PAST:   train past, gate raising
 *
 * @param approach_dist  Distance when train becomes "near"
 * @param crossing_start Start of crossing
 * @param crossing_end   End of crossing
 * @param gate_down      Time to lower gate
 * @param gate_up        Time to raise gate
 * @param train_speed    Train velocity
 * @return               Hybrid automaton
 */
HybridAutomaton* example_train_gate(double approach_dist, double crossing_start,
                                     double crossing_end, double gate_down,
                                     double gate_up, double train_speed)
{
    HybridAutomaton *ha = hybrid_automaton_create("TrainGate", 2);
    if (!ha) return NULL;

    hybrid_automaton_set_variable(ha, 0, "train_x", HAVAR_REAL, -INFINITY, INFINITY);
    hybrid_automaton_set_variable(ha, 1, "gate_angle", HAVAR_REAL_BOUNDED, 0.0, 90.0);

    /* Train position flows rightward at constant speed */
    /* Gate angle: 0 = up, 90 = down */

    /* Mode 0: FAR (train approaching, gate up) */
    int far = hybrid_automaton_add_mode(ha, "FAR", HAMODE_INITIAL);
    double A[4] = {0, 0, 0, 0};
    double b_far[2] = {train_speed, 0};
    hybrid_flow_set(ha, far, HAFLOW_AFFINE, A, b_far);

    /* Invariant: train_x ≤ approach_dist (1 constr × 2 vars = 2 entries) */
    double H_far[2] = {1, 0};
    double k_far = approach_dist;
    hybrid_invariant_set(ha, far, 1, H_far, &k_far);

    /* Mode 1: NEAR (gate lowering) */
    int near = hybrid_automaton_add_mode(ha, "NEAR", HAMODE_NORMAL);
    double b_near[2] = {train_speed, 90.0 / gate_down};
    hybrid_flow_set(ha, near, HAFLOW_AFFINE, A, b_near);

    /* Invariant: approach_dist < train_x < crossing_start, gate_angle < 90
       3 constraints × 2 vars = 6 entries, row-major */
    double H_near[6] = {-1, 0,  1, 0,  0, 1};
    double k_near[3] = {-approach_dist, crossing_start, 90};
    hybrid_invariant_set(ha, near, 3, H_near, k_near);

    /* Mode 2: IN (train crossing, gate down) */
    int in_mode = hybrid_automaton_add_mode(ha, "IN", HAMODE_NORMAL);
    double b_in[2] = {train_speed, 0};
    hybrid_flow_set(ha, in_mode, HAFLOW_AFFINE, A, b_in);

    /* Invariant: crossing_start ≤ train_x ≤ crossing_end
       2 constraints × 2 vars = 4 entries, row-major */
    double H_in[4] = {-1, 0,  1, 0};
    double k_in[2] = {-crossing_start, crossing_end};
    hybrid_invariant_set(ha, in_mode, 2, H_in, k_in);

    /* Mode 3: PAST (gate raising) */
    int past = hybrid_automaton_add_mode(ha, "PAST", HAMODE_NORMAL);
    double b_past[2] = {train_speed, -90.0 / gate_up};
    hybrid_flow_set(ha, past, HAFLOW_AFFINE, A, b_past);

    /* Transitions: FAR→NEAR, NEAR→IN, IN→PAST */
    int t1 = hybrid_automaton_add_transition(ha, far, near, HATRIG_AUTONOMOUS, "approach");
    hybrid_guard_set(ha, t1, 1, H_far, &k_far);
    hybrid_reset_set(ha, t1, HARESET_IDENTITY, NULL, NULL);

    int t2 = hybrid_automaton_add_transition(ha, near, in_mode, HATRIG_AUTONOMOUS, "enter");
    hybrid_guard_set(ha, t2, 1, (double[]){1}, (double[]){crossing_start});
    hybrid_reset_set(ha, t2, HARESET_IDENTITY, NULL, NULL);

    int t3 = hybrid_automaton_add_transition(ha, in_mode, past, HATRIG_AUTONOMOUS, "exit");
    hybrid_guard_set(ha, t3, 1, (double[]){1}, (double[]){crossing_end});
    hybrid_reset_set(ha, t3, HARESET_IDENTITY, NULL, NULL);

    double x0[2] = {0.0, 0.0};
    hybrid_init_set(ha, far, x0);

    return ha;
}

/* ==========================================================================
 * KP53: DC-DC Converter
 * ========================================================================== */

/**
 * @brief KP53: Create a DC-DC boost converter hybrid automaton.
 *
 * Two modes: switch ON (inductor charges) and switch OFF (inductor discharges
 * to capacitor). State includes inductor current i_L and capacitor voltage v_C.
 *
 * State: x = (i_L, v_C)
 *
 * @param L    Inductance
 * @param C    Capacitance
 * @param R    Load resistance
 * @param V_in Input voltage
 * @param iL0  Initial inductor current
 * @param vC0  Initial capacitor voltage
 * @return     Hybrid automaton
 */
HybridAutomaton* example_dcdc_converter(double L, double C, double R,
                                         double V_in, double iL0, double vC0)
{
    HybridAutomaton *ha = hybrid_automaton_create("DCDC_Boost", 2);
    if (!ha) return NULL;

    hybrid_automaton_set_variable(ha, 0, "iL", HAVAR_REAL, 0.0, INFINITY);
    hybrid_automaton_set_variable(ha, 1, "vC", HAVAR_REAL, 0.0, INFINITY);

    /* Mode 0: SWITCH ON
       diL/dt = V_in / L
       dvC/dt = -vC / (R·C) */
    int mode_on = hybrid_automaton_add_mode(ha, "SW_ON", HAMODE_INITIAL);
    double A_on[4] = {0, 0, 0, -1.0 / (R * C)};
    double b_on[2] = {V_in / L, 0};
    hybrid_flow_set(ha, mode_on, HAFLOW_AFFINE, A_on, b_on);

    /* Mode 1: SWITCH OFF
       diL/dt = (V_in - v_C) / L
       dvC/dt = (i_L - v_C/R) / C */
    int mode_off = hybrid_automaton_add_mode(ha, "SW_OFF", HAMODE_NORMAL);
    double A_off[4] = {0, -1.0 / L, 1.0 / C, -1.0 / (R * C)};
    double b_off[2] = {V_in / L, 0};
    hybrid_flow_set(ha, mode_off, HAFLOW_AFFINE, A_off, b_off);

    /* Controlled transitions (switch control) */
    int on_to_off = hybrid_automaton_add_transition(ha, mode_on, mode_off,
                                                      HATRIG_CONTROLLED, "switch_off");
    hybrid_reset_set(ha, on_to_off, HARESET_IDENTITY, NULL, NULL);

    int off_to_on = hybrid_automaton_add_transition(ha, mode_off, mode_on,
                                                      HATRIG_CONTROLLED, "switch_on");
    hybrid_reset_set(ha, off_to_on, HARESET_IDENTITY, NULL, NULL);

    double x0[2] = {iL0, vC0};
    hybrid_init_set(ha, mode_on, x0);

    return ha;
}

/* ==========================================================================
 * KP54: Engine Air-Fuel Ratio Control (L7 Application)
 * ========================================================================== */

/**
 * @brief KP54: Engine air-fuel ratio (AFR) hybrid control model.
 *
 * Automotive application: the AFR cycles between lean and rich
 * around stoichiometric (λ = 1) for optimal catalytic converter
 * efficiency.
 *
 * State: x = (AFR, catalyst_O2_storage)
 *
 * Modes: LEAN (excess air), RICH (excess fuel)
 * Transitions: based on oxygen sensor readings
 *
 * @return Hybrid automaton
 */
HybridAutomaton* example_engine_afr_control(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("EngineAFR", 2);
    if (!ha) return NULL;

    hybrid_automaton_set_variable(ha, 0, "AFR", HAVAR_REAL_BOUNDED, 10.0, 20.0);
    hybrid_automaton_set_variable(ha, 1, "O2_storage", HAVAR_REAL_BOUNDED, 0.0, 1.0);

    /* Mode 0: LEAN (λ > 1) */
    int lean = hybrid_automaton_add_mode(ha, "LEAN", HAMODE_INITIAL);
    double A_lean[4] = {-5, 0, 0, 0.1};
    double b_lean[2] = {75, 0};
    hybrid_flow_set(ha, lean, HAFLOW_AFFINE, A_lean, b_lean);

    /* Mode 1: RICH (λ < 1) */
    int rich = hybrid_automaton_add_mode(ha, "RICH", HAMODE_NORMAL);
    double A_rich[4] = {-5, 0, 0, -0.2};
    double b_rich[2] = {65, 0};
    hybrid_flow_set(ha, rich, HAFLOW_AFFINE, A_rich, b_rich);

    /* Transitions */
    int l2r = hybrid_automaton_add_transition(ha, lean, rich, HATRIG_CONTROLLED, "go_rich");
    int r2l = hybrid_automaton_add_transition(ha, rich, lean, HATRIG_CONTROLLED, "go_lean");
    hybrid_reset_set(ha, l2r, HARESET_IDENTITY, NULL, NULL);
    hybrid_reset_set(ha, r2l, HARESET_IDENTITY, NULL, NULL);

    double x0[2] = {14.7, 0.5};
    hybrid_init_set(ha, lean, x0);

    return ha;
}

/* ==========================================================================
 * KP55: Robot Obstacle Avoidance (L7 Application)
 * ========================================================================== */

/**
 * @brief KP55: Robot navigation with obstacle avoidance.
 *
 * Mobile robot in 2D plane, obstacle detection triggers mode switches
 * between cruising (straight line) and avoidance (circular arc).
 *
 * State: x = (pos_x, pos_y, heading_theta)
 *
 * Modes: CRUISE, AVOID_LEFT, AVOID_RIGHT, REJOIN
 *
 * @param cruise_speed Cruising speed
 * @param turn_radius  Avoidance turn radius
 * @param obstacle_x, obstacle_y Obstacle position
 * @return Hybrid automaton
 */
HybridAutomaton* example_robot_obstacle_avoidance(double cruise_speed,
                                                    double turn_radius,
                                                    double obstacle_x,
                                                    double obstacle_y)
{
    HybridAutomaton *ha = hybrid_automaton_create("RobotAvoid", 3);
    if (!ha) return NULL;

    /* Store obstacle position for reference (application-level data) */
    (void)obstacle_x;
    (void)obstacle_y;

    hybrid_automaton_set_variable(ha, 0, "x", HAVAR_REAL, -INFINITY, INFINITY);
    hybrid_automaton_set_variable(ha, 1, "y", HAVAR_REAL, -INFINITY, INFINITY);
    hybrid_automaton_set_variable(ha, 2, "theta", HAVAR_REAL, -M_PI, M_PI);

    /* Mode 0: CRUISE (straight line) */
    int cruise = hybrid_automaton_add_mode(ha, "CRUISE", HAMODE_INITIAL);
    double b_cruise[3] = {cruise_speed, 0, 0};
    hybrid_flow_set(ha, cruise, HAFLOW_CONSTANT, NULL, b_cruise);

    /* Mode 1: AVOID_LEFT */
    int avoid_l = hybrid_automaton_add_mode(ha, "AVOID_L", HAMODE_NORMAL);
    double omega = cruise_speed / turn_radius;
    /* ẋ = v·cos(θ), ẏ = v·sin(θ), θ̇ = ω */
    /* Nonlinear approximation using constant values at initial θ */
    double b_avoid_l[3] = {cruise_speed * 0.7, cruise_speed * 0.7, omega};
    hybrid_flow_set(ha, avoid_l, HAFLOW_CONSTANT, NULL, b_avoid_l);

    /* Mode 2: REJOIN (return to path) */
    int rejoin = hybrid_automaton_add_mode(ha, "REJOIN", HAMODE_NORMAL);
    double b_rejoin[3] = {cruise_speed, 0, -omega * 0.5};
    hybrid_flow_set(ha, rejoin, HAFLOW_CONSTANT, NULL, b_rejoin);

    /* Transitions */
    int c2a = hybrid_automaton_add_transition(ha, cruise, avoid_l,
                                                HATRIG_CONTROLLED, "detect_obstacle");
    hybrid_reset_set(ha, c2a, HARESET_IDENTITY, NULL, NULL);

    int a2r = hybrid_automaton_add_transition(ha, avoid_l, rejoin,
                                                HATRIG_CONTROLLED, "cleared");
    hybrid_reset_set(ha, a2r, HARESET_IDENTITY, NULL, NULL);

    int r2c = hybrid_automaton_add_transition(ha, rejoin, cruise,
                                                HATRIG_CONTROLLED, "rejoined");
    hybrid_reset_set(ha, r2c, HARESET_IDENTITY, NULL, NULL);

    double x0[3] = {0, 0, 0};
    hybrid_init_set(ha, cruise, x0);

    return ha;
}

/* ==========================================================================
 * KP56: Medical Infusion Pump (L7 Application)
 * ========================================================================== */

/**
 * @brief KP56: Medical infusion pump safety model.
 *
 * Infusion pump delivering medication with flow rate control.
 * Safety requirement: drug concentration stays within therapeutic
 * window.
 *
 * State: x = (drug_concentration, pump_rate)
 *
 * Modes: NORMAL_FLOW, HIGH_ALARM, LOW_ALARM, STOPPED
 *
 * @return Hybrid automaton
 */
HybridAutomaton* example_infusion_pump(void)
{
    HybridAutomaton *ha = hybrid_automaton_create("InfusionPump", 2);
    if (!ha) return NULL;

    hybrid_automaton_set_variable(ha, 0, "conc", HAVAR_REAL_BOUNDED, 0.0, 100.0);
    hybrid_automaton_set_variable(ha, 1, "rate", HAVAR_REAL_BOUNDED, 0.0, 10.0);

    /* Mode 0: NORMAL */
    int normal = hybrid_automaton_add_mode(ha, "NORMAL", HAMODE_INITIAL);
    double A_norm[4] = {-0.1, 1, 0, 0};
    double b_norm[2] = {0, 0};
    hybrid_flow_set(ha, normal, HAFLOW_LINEAR, A_norm, b_norm);

    /* Invariant: 0 ≤ conc ≤ 80 (2 constraints × 2 vars = 4 entries) */
    double H_norm[4] = {1, 0,  -1, 0};
    double k_norm[2] = {80, 0};
    hybrid_invariant_set(ha, normal, 2, H_norm, k_norm);

    /* Mode 1: HIGH_ALARM */
    int high = hybrid_automaton_add_mode(ha, "HIGH_ALARM", HAMODE_NORMAL);
    double b_high[2] = {0, 0};
    hybrid_flow_set(ha, high, HAFLOW_CONSTANT, NULL, b_high);

    /* Mode 2: STOPPED */
    int stopped = hybrid_automaton_add_mode(ha, "STOPPED", HAMODE_TERMINAL);
    double b_stop[2] = {0, 0};
    hybrid_flow_set(ha, stopped, HAFLOW_CONSTANT, NULL, b_stop);

    /* Transitions */
    int n2h = hybrid_automaton_add_transition(ha, normal, high, HATRIG_AUTONOMOUS, "conc_high");
    double gh[1] = {1};
    double gb_h = 80;
    hybrid_guard_set(ha, n2h, 1, gh, &gb_h);
    hybrid_reset_set(ha, n2h, HARESET_CONSTANT, NULL, (double[]){0, 0});

    int h2n = hybrid_automaton_add_transition(ha, high, normal, HATRIG_CONTROLLED, "acknowledge");
    hybrid_reset_set(ha, h2n, HARESET_IDENTITY, NULL, NULL);

    int h2s = hybrid_automaton_add_transition(ha, high, stopped, HATRIG_CONTROLLED, "stop");
    hybrid_reset_set(ha, h2s, HARESET_IDENTITY, NULL, NULL);

    double x0[2] = {0, 5};
    hybrid_init_set(ha, normal, x0);

    return ha;
}
