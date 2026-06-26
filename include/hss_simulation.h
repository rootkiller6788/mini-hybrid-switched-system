/**
 * @file hss_simulation.h
 * @brief Hybrid Switched System Simulation Engine (L5-L6)
 *
 * Provides numerical simulation of hybrid switched systems,
 * implementing event detection, mode switching, and continuous
 * integration within each mode.
 *
 * Knowledge Points:
 *   L5 KP9: Hybrid simulation loop (event detection + ODE stepping)
 *   L5 KP10: Guard crossing detection via root-finding
 *   L5 KP11: Adaptive step-size control with switching
 *   L5 KP12: Multi-rate simulation for stiff systems
 *   L6 KP1: Bouncing ball — canonical hybrid automaton
 *   L6 KP2: Thermostat — discrete switching with hysteresis
 *   L6 KP3: DC-DC converter — power electronics switching
 *   L6 KP4: Vehicle cruise control — automotive hybrid system
 *
 * Course Mapping:
 *   MIT 6.241 — Numerical methods for dynamical systems
 *   Stanford CS205L — Scientific computing
 *   Berkeley EECS 291E — Hybrid system simulation
 *   CMU 15-424 — CPS simulation
 *   ETH 263-0006 — Computer architecture (simulation)
 */

#ifndef HSS_SIMULATION_H
#define HSS_SIMULATION_H

#include "hss_core.h"

/* ============================================================================
 * L5 KP9: Simulation Event
 * ============================================================================ */

/**
 * @brief Simulation event record
 *
 * Records discrete events during hybrid simulation.
 */
typedef struct {
    HSS_EventType  type;          /**< Event type                       */
    double         time;          /**< Time of event t_k                 */
    int            src_mode;      /**< Source mode before event          */
    int            dst_mode;      /**< Destination mode after event      */
    int            jump_index;    /**< Discrete jump counter j           */
    double        *pre_state;     /**< State x(t_k⁻) before reset       */
    double        *post_state;    /**< State x(t_k⁺) after reset        */
    HSS_ResetType  reset_type;    /**< Type of reset applied             */
    char           description[HSS_NAME_LEN]; /**< Event description     */
} HSS_SimEvent;

/**
 * @brief Simulation statistics
 */
typedef struct {
    int     num_events;          /**< Total discrete events              */
    int     num_steps;           /**< Total integration steps            */
    int     num_rejected_steps;  /**< Steps rejected (adaptive)         */
    int     num_mode_switches;   /**< Number of mode changes            */
    double  sim_time;            /**< Total simulation time (wall clock)*/
    double  cpu_time;            /**< CPU time used                     */
    double  min_step;            /**< Minimum step size used            */
    double  max_step;            /**< Maximum step size used            */
    double  avg_step;            /**< Average step size                 */
    double  total_ode_time;      /**< Time spent in ODE solver          */
    double  total_event_time;    /**< Time spent in event detection     */
} HSS_SimStats;

/* ============================================================================
 * L5 KP10-KP12: Advanced Simulation Features
 * ============================================================================ */

/**
 * @brief Root-finding context for guard crossing detection (L5 KP10)
 *
 * Uses zero-crossing detection: find t* where g(x(t*)) = 0
 * for guard function g(x) = aᵀx - b.
 * Implements Illinois algorithm (modified regula falsi).
 */
typedef struct {
    double  t_left;           /**< Left bracket time                   */
    double  t_right;          /**< Right bracket time                  */
    double  g_left;           /**< Guard value at left                 */
    double  g_right;          /**< Guard value at right                */
    int     iterations;       /**< Root-finding iterations              */
    double  t_root;           /**< Found root (event time)             */
    double  g_root;           /**< Guard value at root (≈0)            */
    bool    converged;        /**< True if root found                   */
    int     affected_guard;   /**< Index of crossing guard             */
} HSS_RootFinder;

/**
 * @brief Adaptive step controller (L5 KP11)
 *
 * PI-controller for step size adaptation:
 *   h_{new} = h * (tol / err)^{k_I · k_P} * (err_{old} / err)^{k_P}
 * where k_P = 0.4/p, k_I = 0.3/p for method of order p.
 */
typedef struct {
    double  h_current;        /**< Current step size                    */
    double  h_min;            /**< Minimum allowed step size            */
    double  h_max;            /**< Maximum allowed step size            */
    double  abstol;           /**< Absolute tolerance                   */
    double  reltol;           /**< Relative tolerance                   */
    double  err_prev;         /**< Previous error estimate              */
    double  safety_factor;    /**< Safety factor (typically 0.9)        */
    double  kP;               /**< PI proportional gain                 */
    double  kI;               /**< PI integral gain                     */
    int     order;            /**< Method order p                       */
    int     n_accepted;       /**< Accepted steps                       */
    int     n_rejected;       /**< Rejected steps                       */
} HSS_StepController;

/**
 * @brief Multi-rate simulation context (L5 KP12)
 *
 * For stiff hybrid systems: fast continuous dynamics in some modes
 * coupled with slow discrete switching. Uses different step sizes
 * for different subsystems.
 */
typedef struct {
    HSS_System  *fast_sys;      /**< Fast dynamics subsystem             */
    HSS_System  *slow_sys;      /**< Slow dynamics subsystem             */
    double       fast_h;        /**< Step size for fast dynamics         */
    double       slow_h;        /**< Step size for slow dynamics         */
    int          fast_steps_per_slow; /**< Fast steps per one slow step  */
    double      *coupling_input;/**< Input from slow to fast             */
    double      *coupling_output;/**< Output from fast to slow           */
    bool         is_synchronized;/**< Synchronization flag               */
} HSS_MultiRateSim;

/* ============================================================================
 * L6 KP1-KP4: Canonical Hybrid Systems
 * ============================================================================ */

/**
 * @brief Bouncing ball model (L6 KP1)
 *
 * Classic hybrid automaton from literature:
 *   Mode: free fall — ẍ = -g (gravity)
 *   Guard: x ≤ 0 (ground impact)
 *   Reset: ẋ⁺ = -c · ẋ⁻ (restitution coefficient 0 < c < 1)
 *
 * Demonstrates: Zeno behavior, event detection, reset maps.
 * This is the "hello world" of hybrid systems.
 */
typedef struct {
    double   gravity;          /**< Gravitational acceleration g (m/s²)  */
    double   restitution;      /**< Coefficient of restitution c ∈ (0,1) */
    double   height;           /**< Current height x (m)                 */
    double   velocity;         /**< Current velocity ẋ (m/s)            */
    double   time;             /**< Elapsed time (s)                    */
    int      bounce_count;     /**< Number of bounces                   */
    double   energy;           /**< Mechanical energy (J/kg)             */
    bool     at_rest;          /**< True if ball has stopped             */
} HSS_BouncingBall;

/**
 * @brief Thermostat with hysteresis (L6 KP2)
 *
 * Two-mode hybrid system:
 *   Mode ON (heating): Ṫ = α(TH - T)
 *   Mode OFF (cooling): Ṫ = -β(T - TL)
 *   Switch ON→OFF when T ≥ T_high
 *   Switch OFF→ON when T ≤ T_low
 *
 * Demonstrates: hysteresis switching, invariant sets, deadband.
 */
typedef struct {
    double   temperature;      /**< Current temperature T                */
    double   T_low;            /**< Lower threshold (switch ON)          */
    double   T_high;           /**< Upper threshold (switch OFF)         */
    double   alpha;            /**< Heating rate coefficient             */
    double   beta;             /**< Cooling rate coefficient             */
    double   T_ambient_hot;    /**< Ambient temp when heating (TH)       */
    double   T_ambient_cold;   /**< Ambient temp when cooling (TL)      */
    bool     heater_on;        /**< True if in ON mode                   */
    double   time;             /**< Elapsed time                         */
    double   duty_cycle;       /**< Fraction of time in ON mode         */
    double   avg_temperature;  /**< Time-averaged temperature            */
} HSS_Thermostat;

/**
 * @brief DC-DC Boost Converter (L6 KP3)
 *
 * Two-mode switched affine system:
 *   Mode 1 (switch ON):  L · diL/dt = Vin, C · dvC/dt = -vC/R
 *   Mode 2 (switch OFF): L · diL/dt = Vin - vC, C · dvC/dt = iL - vC/R
 *
 * Demonstrates: PWM switching, limit cycle, average model.
 * State x = [iL, vC] (inductor current, capacitor voltage).
 */
typedef struct {
    double   iL;               /**< Inductor current (A)                 */
    double   vC;               /**< Capacitor voltage (V)                */
    double   Vin;              /**< Input voltage (V)                    */
    double   L;                /**< Inductance (H)                       */
    double   C;                /**< Capacitance (F)                      */
    double   R;                /**< Load resistance (Ω)                  */
    double   duty_cycle;       /**< PWM duty cycle D ∈ [0,1]             */
    double   freq;             /**< Switching frequency (Hz)             */
    double   time;             /**< Elapsed time                         */
    int      switch_state;     /**< 0=OFF, 1=ON                         */
    double   vref;             /**< Reference output voltage             */
    double   output_ripple;    /**< Peak-to-peak ripple voltage          */
} HSS_DCDCConverter;

/**
 * @brief Adaptive Cruise Control (L6 KP4)
 *
 * Hybrid system with mode-dependent dynamics:
 *   Mode CRUISE: maintain set speed v_set
 *   Mode FOLLOW: maintain safe distance from lead vehicle
 *   Mode BRAKE: emergency deceleration
 *
 * State: x = [position, velocity] for ego vehicle.
 * Demonstrates: state-dependent switching, safety constraints.
 */
typedef struct {
    double   ego_position;     /**< Ego vehicle position (m)             */
    double   ego_velocity;     /**< Ego vehicle speed (m/s)              */
    double   lead_position;    /**< Lead vehicle position (m)            */
    double   lead_velocity;    /**< Lead vehicle speed (m/s)             */
    double   v_set;            /**< Desired cruise speed (m/s)           */
    double   time_headway;     /**< Safe time headway (s)                */
    double   min_gap;          /**< Minimum gap at standstill (m)        */
    double   max_accel;        /**< Maximum acceleration (m/s²)          */
    double   max_decel;        /**< Maximum deceleration (m/s²)          */
    double   time;             /**< Elapsed time                         */
    int      mode;             /**< 0=CRUISE, 1=FOLLOW, 2=BRAKE         */
    double   gap_error;        /**< Actual gap - desired gap             */
    bool     collision_risk;   /**< True if collision risk detected      */
} HSS_CruiseControl;

/* ============================================================================
 * L5-L6 API: Simulation Functions
 * ============================================================================ */

/* ---- Initialization ---- */

/**
 * @brief Create default simulation configuration
 * @return Default HSS_SimConfig
 */
HSS_SimConfig hss_sim_config_default(void);

/**
 * @brief Initialize step controller (L5 KP11)
 * @param order Method order (4 for RK4)
 * @param abstol Absolute tolerance
 * @param reltol Relative tolerance
 * @return Initialized step controller
 */
HSS_StepController hss_step_controller_init(int order,
                                             double abstol, double reltol);

/* ---- Core simulation ---- */

/**
 * @brief Simulate a hybrid switched system (L5 KP9)
 *
 * Main hybrid simulation loop:
 * 1. While t < t_end:
 *    a. Integrate continuous dynamics in current mode
 *    b. Check for event triggers (guard crossing, invariant violation)
 *    c. If event detected: apply reset, switch mode, increment jump counter
 *    d. Check for Zeno (too many jumps in short time)
 * 2. Record execution trace
 *
 * @param sys HSS system (initialized with initial state)
 * @param config Simulation configuration
 * @param trace Output execution trace (caller allocates)
 * @param stats Output simulation statistics
 * @return 0 on success, -1 on error, -2 on Zeno detected
 *
 * Complexity: O(n³·steps + events·n) for linear modes,
 *             O(eval(flow)·steps) for nonlinear
 */
int hss_simulate(HSS_System *sys,
                  const HSS_SimConfig *config,
                  HSS_ExecutionTrace *trace,
                  HSS_SimStats *stats);

/**
 * @brief Single step of continuous integration
 *
 * Advances continuous state from t to t+h within current mode,
 * using the configured solver.
 *
 * @param sys HSS system
 * @param h Step size
 * @return 0 on success, -1 on error
 */
int hss_step_continuous(HSS_System *sys, double h);

/**
 * @brief Check all guards for the current mode (L5 KP10)
 *
 * Evaluates all outgoing transition guards and returns the
 * index of the first satisfied transition.
 *
 * @param sys HSS system
 * @return Transition index if guard satisfied, -1 if none
 */
int hss_check_guards(const HSS_System *sys);

/**
 * @brief Apply a discrete transition
 *
 * Executes a specific transition: validates guard, applies reset map,
 * switches mode, increments jump counter.
 *
 * @param sys HSS system
 * @param transition_id Index of transition to fire
 * @return 0 on success, -1 if transition not valid
 */
int hss_apply_transition(HSS_System *sys, int transition_id);

/**
 * @brief Find guard crossing time via root-finding (L5 KP10)
 *
 * Uses Illinois algorithm to find exact crossing time of a guard.
 *
 * @param sys HSS system
 * @param guard_index Index of guard to check
 * @param t_now Current time
 * @param t_next Proposed next time
 * @param rf Output root-finder context
 * @return 0 on success, -1 if no crossing found
 */
int hss_find_guard_crossing(const HSS_System *sys, int guard_index,
                             double t_now, double t_next,
                             HSS_RootFinder *rf);

/* ---- Adaptive step-size control (L5 KP11) ---- */

/**
 * @brief Adapt step size based on error estimate
 *
 * @param sc Step controller
 * @param err Current error estimate
 * @return New recommended step size
 */
double hss_adapt_step(HSS_StepController *sc, double err);

/**
 * @brief Check if step needs to be shortened to hit event
 * @param h Proposed step
 * @param t_next_event Time of next known event
 * @return Adjusted step size (≤ h)
 */
double hss_adjust_step_for_event(double h, double t_next_event);

/* ---- Multi-rate simulation (L5 KP12) ---- */

/**
 * @brief Initialize multi-rate simulation
 * @param fast Fast dynamics subsystem
 * @param slow Slow dynamics subsystem
 * @param ratio Fast steps per slow step
 * @return Initialized multi-rate context
 */
HSS_MultiRateSim hss_multirate_init(HSS_System *fast, HSS_System *slow,
                                     int ratio);

/**
 * @brief Perform one macro-step of multi-rate simulation
 * @param mrs Multi-rate simulation context
 * @return 0 on success
 */
int hss_multirate_macro_step(HSS_MultiRateSim *mrs);

/* ---- L6 Canonical Problem Solvers ---- */

/**
 * @brief Initialize bouncing ball model (L6 KP1)
 * @param gravity Gravity g (m/s²)
 * @param restitution Restitution coefficient c ∈ (0,1)
 * @param initial_height Starting height (m)
 * @param initial_velocity Starting velocity (m/s)
 * @return Initialized bouncing ball
 */
HSS_BouncingBall hss_bouncing_ball_init(double gravity, double restitution,
                                         double initial_height,
                                         double initial_velocity);

/**
 * @brief Simulate bouncing ball for a given time
 * @param bb Ball state (updated in place)
 * @param dt Time step
 * @return Number of bounces in this step
 */
int hss_bouncing_ball_step(HSS_BouncingBall *bb, double dt);

/**
 * @brief Simulate bouncing ball until it comes to rest
 * @param bb Ball state
 * @param max_time Maximum simulation time
 * @return Total bounce count
 */
int hss_bouncing_ball_simulate(HSS_BouncingBall *bb, double max_time);

/**
 * @brief Initialize thermostat model (L6 KP2)
 * @param T_start Initial temperature
 * @param T_low Low threshold
 * @param T_high High threshold
 * @param alpha Heating coefficient
 * @param beta Cooling coefficient
 * @param T_hot Ambient for heating
 * @param T_cold Ambient for cooling
 * @return Initialized thermostat
 */
HSS_Thermostat hss_thermostat_init(double T_start, double T_low,
                                    double T_high, double alpha,
                                    double beta, double T_hot, double T_cold);

/**
 * @brief Advance thermostat simulation
 * @param therm Thermostat state (updated)
 * @param dt Time step
 * @return 0 if no switch, 1 if switched ON→OFF, -1 if switched OFF→ON
 */
int hss_thermostat_step(HSS_Thermostat *therm, double dt);

/**
 * @brief Run thermostat simulation
 * @param therm Thermostat state
 * @param sim_time Total simulation time
 * @param dt Time step
 */
void hss_thermostat_simulate(HSS_Thermostat *therm,
                              double sim_time, double dt);

/**
 * @brief Initialize DC-DC converter model (L6 KP3)
 * @param Vin Input voltage
 * @param L Inductance
 * @param C Capacitance
 * @param R Load resistance
 * @param freq Switching frequency
 * @param duty_cycle PWM duty cycle
 * @return Initialized converter
 */
HSS_DCDCConverter hss_dcdc_init(double Vin, double L, double C,
                                 double R, double freq, double duty_cycle);

/**
 * @brief Advance DC-DC converter simulation
 * @param conv Converter state (updated)
 * @param dt Time step
 * @return 0 if no switch, 1 if switched
 */
int hss_dcdc_step(HSS_DCDCConverter *conv, double dt);

/**
 * @brief Run DC-DC converter to steady state
 * @param conv Converter state
 * @param sim_time Simulation time
 * @param dt Step size
 * @return Steady-state output voltage
 */
double hss_dcdc_simulate(HSS_DCDCConverter *conv,
                          double sim_time, double dt);

/**
 * @brief Initialize cruise control model (L6 KP4)
 * @param ego_vel Initial velocity
 * @param lead_pos Lead vehicle position
 * @param lead_vel Lead vehicle velocity
 * @param v_set Desired cruise speed
 * @param headway Time headway
 * @param min_gap Minimum gap
 * @return Initialized cruise control
 */
HSS_CruiseControl hss_cruise_control_init(double ego_vel,
                                           double lead_pos, double lead_vel,
                                           double v_set, double headway,
                                           double min_gap);

/**
 * @brief Advance cruise control simulation
 * @param cc Cruise control state (updated)
 * @param dt Time step
 * @return Current mode (0=CRUISE, 1=FOLLOW, 2=BRAKE)
 */
int hss_cruise_control_step(HSS_CruiseControl *cc, double dt);

/**
 * @brief Run cruise control simulation
 * @param cc Cruise control state
 * @param sim_time Simulation time
 * @param dt Step size
 */
void hss_cruise_control_simulate(HSS_CruiseControl *cc,
                                  double sim_time, double dt);

/* ---- Execution trace management ---- */

/**
 * @brief Allocate execution trace
 * @param max_steps Maximum number of time points
 * @param state_dim State dimension
 * @return Allocated trace (caller frees via hss_trace_free)
 */
HSS_ExecutionTrace *hss_trace_alloc(int max_steps, int state_dim);

/**
 * @brief Free execution trace
 * @param trace Trace to free
 */
void hss_trace_free(HSS_ExecutionTrace *trace);

/**
 * @brief Record a data point in the trace
 * @param trace Execution trace
 * @param t Current time
 * @param mode Current mode
 * @param x Current state vector
 * @param n State dimension
 * @return 0 on success, -1 if trace is full
 */
int hss_trace_record(HSS_ExecutionTrace *trace, double t,
                      int mode, const double *x, int n);

/**
 * @brief Print execution trace summary
 * @param trace Execution trace
 * @param fp Output stream
 */
void hss_trace_print_summary(const HSS_ExecutionTrace *trace, FILE *fp);

/**
 * @brief Export trace data to CSV file
 * @param trace Execution trace
 * @param filename Output file path
 * @return 0 on success, -1 on I/O error
 */
int hss_trace_export_csv(const HSS_ExecutionTrace *trace,
                          const char *filename);

#endif /* HSS_SIMULATION_H */
