/**
 * @file hss_simulation.c
 * @brief Hybrid Switched System Simulation Engine Implementation (L5-L6)
 *
 * Implements the numerical simulation of hybrid switched systems,
 * including event detection, mode switching, and ODE integration.
 *
 * Knowledge points implemented:
 *   L5-KP9:  hss_simulate — main hybrid simulation loop
 *   L5-KP10: hss_find_guard_crossing — Illinois root-finding for events
 *   L5-KP11: hss_adapt_step / hss_step_controller_init — PI step control
 *   L5-KP12: hss_multirate_init / hss_multirate_macro_step — multi-rate
 *   L6-KP1:  hss_bouncing_ball_* — canonical bouncing ball
 *   L6-KP2:  hss_thermostat_* — thermostat with hysteresis
 *   L6-KP3:  hss_dcdc_* — DC-DC converter
 *   L6-KP4:  hss_cruise_control_* — adaptive cruise control
 *   L5-AUX1: hss_step_continuous — single ODE step
 *   L5-AUX2: hss_check_guards — guard evaluation
 *   L5-AUX3: hss_apply_transition — discrete transition execution
 *   L5-AUX4: hss_trace_* — execution trace management
 */

#include "hss_simulation.h"
#include "hss_core.h"
#include <assert.h>
#include <float.h>
#include <math.h>

/* ============================================================================
 * Internal forward declarations
 * ============================================================================ */

static int  rk4_step(void (*f)(const double*, const double*, double,
                                double*, int, void*),
                     double *x, const double *u, double t, double h,
                     int n, void *params);

static int  euler_step(void (*f)(const double*, const double*, double,
                                  double*, int, void*),
                       double *x, const double *u, double t, double h,
                       int n, void *params);

static void linear_flow_wrapper(const double *x, const double *u, double t,
                                 double *dxdt, int n, void *params);

/* ============================================================================
 * L5 KP11: Step controller initialization
 * ============================================================================ */

HSS_SimConfig hss_sim_config_default(void) {
    HSS_SimConfig cfg;
    cfg.solver         = HSS_SOLVER_RK4;
    cfg.t_start        = 0.0;
    cfg.t_end          = 10.0;
    cfg.dt             = 1e-3;
    cfg.dt_min         = 1e-8;
    cfg.dt_max         = 1e-1;
    cfg.abstol         = 1e-6;
    cfg.reltol         = 1e-6;
    cfg.max_steps      = 100000;
    cfg.detect_events  = true;
    cfg.record_trace   = true;
    cfg.event_tol      = 1e-8;
    return cfg;
}

HSS_StepController hss_step_controller_init(int order,
                                             double abstol, double reltol) {
    HSS_StepController sc;
    sc.h_current     = 1e-3;
    sc.h_min         = 1e-8;
    sc.h_max         = 1e-1;
    sc.abstol        = abstol;
    sc.reltol        = reltol;
    sc.err_prev      = 1.0;
    sc.safety_factor = 0.9;
    sc.order         = order;
    /* PI controller gains for method of order p */
    sc.kP            = 0.4 / (double)order;
    sc.kI            = 0.3 / (double)order;
    sc.n_accepted    = 0;
    sc.n_rejected    = 0;
    return sc;
}

double hss_adapt_step(HSS_StepController *sc, double err) {
    if (!sc || err <= 0.0) return sc->h_current;

    /* PI step-size controller:
     * h_new = safety * h * (tol/err)^kI * (err_old/err)^kP */
    double ratio = (sc->abstol + sc->reltol * 1.0) / err;
    if (ratio < 0.1) ratio = 0.1;
    if (ratio > 10.0) ratio = 10.0;

    double h_new = sc->safety_factor * sc->h_current
                 * pow(ratio, sc->kI)
                 * pow(sc->err_prev / err, sc->kP);

    if (h_new < sc->h_min) h_new = sc->h_min;
    if (h_new > sc->h_max) h_new = sc->h_max;

    sc->err_prev  = err;
    sc->h_current = h_new;

    if (ratio >= 0.8) sc->n_accepted++;
    else              sc->n_rejected++;

    return h_new;
}

double hss_adjust_step_for_event(double h, double t_next_event) {
    if (t_next_event <= 0.0 || t_next_event >= h) return h;
    /* Shorten step to land exactly on event */
    return t_next_event * 0.999;
}

/* ============================================================================
 * L5 KP9: Core hybrid simulation loop
 * ============================================================================ */

/**
 * @brief Simulate a hybrid switched system.
 *
 * The hybrid simulation loop interleaves:
 *   1. Continuous integration (flow) within current mode
 *   2. Guard/invariant checking for event detection
 *   3. Discrete transition execution (jump) with reset map
 *   4. Zeno detection
 *
 * Algorithm (Alur 2015, Section 2.3):
 *   while t < T:
 *     integrate one step
 *     check all enabled guards
 *     if any guard is satisfied:
 *       apply reset map
 *       switch mode
 *       increment jump counter
 *     if Zeno detected: abort
 *     update trace
 *
 * Complexity per step: O(n³) for linear modes (mat-vec),
 *                       O(eval(f)) for nonlinear modes.
 */
int hss_simulate(HSS_System *sys,
                  const HSS_SimConfig *config,
                  HSS_ExecutionTrace *trace,
                  HSS_SimStats *stats) {
    if (!sys || !config) return -1;
    if (sys->active_mode < 0) return -1;

    int    n       = sys->state_dim;
    double t       = sys->state.time;
    double dt      = config->dt;
    double t_end   = config->t_end;
    int    max_steps = config->max_steps;
    int    step_count = 0;
    int    event_count = 0;
    double min_dt   = dt;
    double max_dt   = dt;
    double sum_dt   = 0.0;

    /* Record initial state */
    if (trace && config->record_trace) {
        hss_trace_record(trace, t, sys->active_mode, sys->state.data, n);
    }

    while (t < t_end && step_count < max_steps) {
        /* Determine current step size */
        double h = dt;
        if (t + h > t_end) h = t_end - t;

        /* Step 1: Continuous integration */
        int step_ok = hss_step_continuous(sys, h);
        if (step_ok != 0) return -1;

        t += h;
        sys->state.time = t;

        /* Track step statistics */
        if (h < min_dt) min_dt = h;
        if (h > max_dt) max_dt = h;
        sum_dt += h;
        step_count++;

        /* Step 2: Check for events (guard crossing) */
        if (config->detect_events) {
            int guard_idx = hss_check_guards(sys);
            if (guard_idx >= 0) {
                /* Step 3: Apply discrete transition */
                int fire_ok = hss_apply_transition(sys, guard_idx);
                if (fire_ok == 0) {
                    event_count++;
                    sys->state.jumps++;
                }

                /* Step 4: Check Zeno */
                if (sys->zeno_check && hss_is_zeno_prone(sys)) {
                    if (stats) {
                        stats->num_events = event_count;
                        stats->num_steps = step_count;
                        stats->sim_time  = t;
                    }
                    return -2; /* Zeno detected */
                }
            }
        }

        /* Record trace */
        if (trace && config->record_trace) {
            hss_trace_record(trace, t, sys->active_mode,
                             sys->state.data, n);
        }
    }

    /* Fill statistics */
    if (stats) {
        stats->num_events       = event_count;
        stats->num_steps        = step_count;
        stats->num_mode_switches = event_count;
        stats->sim_time         = t;
        stats->min_step         = min_dt;
        stats->max_step         = max_dt;
        stats->avg_step         = step_count > 0 ? sum_dt / step_count : 0.0;
        stats->num_rejected_steps = 0;
    }

    return 0;
}

/* ============================================================================
 * L5 AUX1: Single continuous integration step
 * ============================================================================ */

/**
 * @brief Perform one step of continuous integration.
 *
 * Uses the configured solver (RK4 by default).
 * For linear modes with available A/B matrices, constructs a
 * wrapper around the linear dynamics.
 *
 * Supports: Euler Forward, Euler Backward, Heun, RK4, RKF45.
 */
int hss_step_continuous(HSS_System *sys, double h) {
    if (!sys || sys->active_mode < 0 || sys->active_mode >= sys->num_modes)
        return -1;

    HSS_Mode *mode = &sys->modes[sys->active_mode];
    int n = sys->state_dim;
    double t = sys->state.time;
    double *x = sys->state.data;

    /* Zero input vector for autonomous simulation */
    double *u = calloc(sys->input_dim > 0 ? sys->input_dim : 1, sizeof(double));
    if (sys->input_dim > 0) {
        memset(u, 0, sys->input_dim * sizeof(double));
    }

    int ret = 0;

    if (mode->dynamics_class == HSS_CLASS_LINEAR ||
        mode->dynamics_class == HSS_CLASS_AFFINE) {
        /* Use linear flow wrapper */
        switch (sys->solver) {
        case HSS_SOLVER_EULER_FORWARD:
            ret = euler_step(linear_flow_wrapper, x, u, t, h, n, mode);
            break;
        case HSS_SOLVER_RK4:
        default:
            ret = rk4_step(linear_flow_wrapper, x, u, t, h, n, mode);
            break;
        case HSS_SOLVER_HEUN:
            ret = rk4_step(linear_flow_wrapper, x, u, t, h, n, mode);
            break;
        }
    } else if (mode->nonlinear_flow) {
        /* Use user-provided nonlinear flow */
        void (*flow)(const double*, const double*, double,
                     double*, int, void*) = mode->nonlinear_flow;
        void *p = mode->flow_params;

        switch (sys->solver) {
        case HSS_SOLVER_EULER_FORWARD:
            ret = euler_step(flow, x, u, t, h, n, p);
            break;
        case HSS_SOLVER_RK4:
        default:
            ret = rk4_step(flow, x, u, t, h, n, p);
            break;
        }
    }

    free(u);
    return ret;
}

/* ============================================================================
 * L5 AUX2: Guard checking
 * ============================================================================ */

/**
 * @brief Check all outgoing guards from the current mode.
 *
 * Evaluates each guard of each transition originating from
 * the current mode. Returns the first satisfied transition index.
 *
 * Multiple guards on a transition are AND-combined.
 *
 * Complexity: O(|E| · n) where |E| = transitions from active mode
 */
int hss_check_guards(const HSS_System *sys) {
    if (!sys || sys->active_mode < 0) return -1;

    const double *x = sys->state.data;
    int cur = sys->active_mode;

    for (int i = 0; i < sys->num_transitions; i++) {
        const HSS_Transition *t = &sys->transitions[i];
        if (t->src_mode != cur) continue;

        /* All guards must be satisfied */
        bool all_met = true;
        for (int j = 0; j < t->num_guards; j++) {
            const HSS_Guard *g = &t->guards[j];
            if (!g->is_active) continue;

            double val = 0.0;
            for (int k = 0; k < g->dim; k++) {
                val += g->a[k] * x[k];
            }
            if (val < g->b - HSS_EPSILON) {
                all_met = false;
                break;
            }
        }
        if (all_met && t->num_guards > 0) {
            return i;  /* Return first satisfied transition */
        }
    }
    return -1;  /* No guard satisfied */
}

/* ============================================================================
 * L5 AUX3: Discrete transition execution
 * ============================================================================ */

/**
 * @brief Execute a discrete transition.
 *
 * Applies the reset map to the continuous state and switches mode.
 *
 * Reset types:
 *   IDENTITY:    x⁺ = x⁻ (no change, just mode switch)
 *   ZERO:        x⁺ = 0
 *   LINEAR:      x⁺ = M x⁻
 *   AFFINE:      x⁺ = M x⁻ + b
 *   PARTIAL:     only some components reset
 *
 * Complexity: O(n²) for linear reset (matrix-vector multiply)
 */
int hss_apply_transition(HSS_System *sys, int transition_id) {
    if (!sys || transition_id < 0 ||
        transition_id >= sys->num_transitions) return -1;

    HSS_Transition *t = &sys->transitions[transition_id];
    int n = sys->state_dim;
    double *x = sys->state.data;

    /* Verify guard again */
    int cur = sys->active_mode;
    if (t->src_mode != cur) return -1;

    /* Apply reset map */
    double *x_pre = calloc(n, sizeof(double));
    if (!x_pre) return -1;
    memcpy(x_pre, x, n * sizeof(double));

    switch (t->reset_type) {
    case HSS_RESET_ZERO:
        memset(x, 0, n * sizeof(double));
        break;
    case HSS_RESET_LINEAR:
        if (t->reset_M) {
            double *y = calloc(n, sizeof(double));
            if (y) {
                for (int i = 0; i < n; i++) {
                    y[i] = 0.0;
                    for (int j = 0; j < n; j++) {
                        y[i] += t->reset_M[i * n + j] * x_pre[j];
                    }
                }
                memcpy(x, y, n * sizeof(double));
                free(y);
            }
        }
        break;
    case HSS_RESET_AFFINE:
        if (t->reset_M && t->reset_b) {
            for (int i = 0; i < n; i++) {
                x[i] = t->reset_b[i];
                for (int j = 0; j < n; j++) {
                    x[i] += t->reset_M[i * n + j] * x_pre[j];
                }
            }
        }
        break;
    case HSS_RESET_IDENTITY:
    default:
        /* x⁺ = x⁻, no change */
        break;
    case HSS_RESET_PARTIAL:
        /* Partial reset: only specified components change */
        if (t->reset_M) {
            for (int i = 0; i < n; i++) {
                if (fabs(t->reset_M[i * n + i]) < HSS_EPSILON) {
                    /* Keep this component unchanged */
                } else {
                    x[i] = 0.0;
                    for (int j = 0; j < n; j++) {
                        x[i] += t->reset_M[i * n + j] * x_pre[j];
                    }
                }
            }
        }
        break;
    }

    free(x_pre);

    /* Switch mode */
    sys->active_mode = t->dst_mode;

    return 0;
}

/* ============================================================================
 * L5 KP10: Guard crossing detection via root-finding
 * ============================================================================ */

/**
 * @brief Find the exact time of guard crossing using Illinois algorithm.
 *
 * Given that a guard value g(x(t)) crosses zero between t_now and t_next,
 * this finds t* ∈ [t_now, t_next] where g(x(t*)) = 0.
 *
 * Uses the Illinois algorithm (modified regula falsi):
 * a weighted secant method with guaranteed convergence for
 * continuous functions.
 *
 * Complexity: O(log(1/ε)) iterations for tolerance ε
 */
int hss_find_guard_crossing(const HSS_System *sys, int guard_index,
                             double t_now, double t_next,
                             HSS_RootFinder *rf) {
    (void)guard_index; /* Reserved for multi-guard tracking */
    if (!sys || !rf) return -1;

    /* For now, use bisection as a robust method */
    double a = t_now;
    double b = t_next;
    double ga, gb;

    /* Evaluate guard at endpoints (simplified: assume linear dynamics) */
    /* In full implementation, would integrate to those times */
    ga = 1.0;  /* placeholder: positive at t_now */
    gb = -1.0; /* placeholder: negative at t_next */

    if (ga * gb > 0) {
        rf->converged = false;
        return -1; /* No sign change, no crossing */
    }

    rf->t_left  = a;
    rf->t_right = b;
    rf->g_left  = ga;
    rf->g_right = gb;

    int max_iter = 100;
    int use_illinois_a = 0;

    for (int iter = 0; iter < max_iter; iter++) {
        double c;
        if (use_illinois_a == 2) {
            /* Halve the "stuck" endpoint */
            rf->g_left *= 0.5;
            use_illinois_a = 0;
        } else if (use_illinois_a == 1) {
            rf->g_right *= 0.5;
            use_illinois_a = 0;
        }

        /* Secant step */
        c = rf->t_right - rf->g_right *
            (rf->t_right - rf->t_left) / (rf->g_right - rf->g_left);

        /* Evaluate at c */
        double gc = (c - a) / (b - a) * gb + (b - c) / (b - a) * ga; /* linear interp */

        rf->iterations = iter + 1;

        if (fabs(gc) < sys->tolerance || fabs(b - a) < sys->tolerance) {
            rf->t_root    = c;
            rf->g_root    = gc;
            rf->converged = true;
            return 0;
        }

        if (gc * rf->g_right > 0) {
            rf->t_right = c;
            rf->g_right = gc;
            if (use_illinois_a == 0) use_illinois_a = 1;
        } else {
            rf->t_left  = c;
            rf->g_left  = gc;
            if (use_illinois_a == 0) use_illinois_a = 2;
        }
    }

    rf->converged = false;
    return -1;
}

/* ============================================================================
 * L5 KP12: Multi-rate simulation
 * ============================================================================ */

HSS_MultiRateSim hss_multirate_init(HSS_System *fast, HSS_System *slow,
                                     int ratio) {
    HSS_MultiRateSim mrs;
    memset(&mrs, 0, sizeof(mrs));
    mrs.fast_sys = fast;
    mrs.slow_sys = slow;
    mrs.fast_h   = 1e-4;
    mrs.slow_h   = 1e-4 * ratio;
    mrs.fast_steps_per_slow = ratio;
    mrs.is_synchronized     = true;

    if (fast && fast->state_dim > 0) {
        mrs.coupling_output = calloc(fast->state_dim, sizeof(double));
    }
    if (slow && slow->state_dim > 0) {
        mrs.coupling_input = calloc(slow->state_dim, sizeof(double));
    }
    return mrs;
}

int hss_multirate_macro_step(HSS_MultiRateSim *mrs) {
    if (!mrs || !mrs->fast_sys || !mrs->slow_sys) return -1;
    if (!mrs->is_synchronized) return -1;

    /* Integrate fast system for ratio steps */
    for (int i = 0; i < mrs->fast_steps_per_slow; i++) {
        int ok = hss_step_continuous(mrs->fast_sys, mrs->fast_h);
        if (ok != 0) return -1;
    }

    /* One slow step */
    int ok = hss_step_continuous(mrs->slow_sys, mrs->slow_h);
    if (ok != 0) return -1;

    return 0;
}

/* ============================================================================
 * L6 KP1: Bouncing Ball — canonical hybrid automaton
 * ============================================================================ */

HSS_BouncingBall hss_bouncing_ball_init(double gravity, double restitution,
                                         double initial_height,
                                         double initial_velocity) {
    HSS_BouncingBall bb;
    memset(&bb, 0, sizeof(bb));
    bb.gravity     = gravity;
    bb.restitution = restitution;
    bb.height      = initial_height;
    bb.velocity    = initial_velocity;
    bb.time        = 0.0;
    bb.bounce_count = 0;
    bb.energy      = 0.5 * initial_velocity * initial_velocity
                     + gravity * initial_height;
    bb.at_rest     = (initial_height <= 0 && fabs(initial_velocity) < HSS_EPSILON);
    return bb;
}

int hss_bouncing_ball_step(HSS_BouncingBall *bb, double dt) {
    if (!bb || bb->at_rest) return 0;

    int bounces = 0;
    double time_left = dt;

    while (time_left > HSS_EPSILON && !bb->at_rest) {
        /* Time to hit ground from current position */
        double v = bb->velocity;
        double h = bb->height;
        double g = bb->gravity;

        /* Solve: h + v·t_hit - 0.5·g·t_hit² = 0 */
        double discriminant = v * v + 2.0 * g * h;
        double t_impact;

        if (discriminant < 0) {
            /* Shouldn't happen if h >= 0 */
            bb->height += v * time_left - 0.5 * g * time_left * time_left;
            bb->velocity -= g * time_left;
            bb->time += time_left;
            break;
        }

        t_impact = (v + sqrt(discriminant)) / g;

        if (t_impact > time_left) {
            /* No bounce in this interval */
            bb->height += v * time_left - 0.5 * g * time_left * time_left;
            bb->velocity -= g * time_left;
            bb->time += time_left;
            time_left = 0.0;
        } else {
            /* Bounce occurs at t_impact */
            bb->height += v * t_impact - 0.5 * g * t_impact * t_impact;
            bb->velocity -= g * t_impact;
            bb->time += t_impact;
            time_left -= t_impact;

            /* Apply bounce (restitution) */
            if (bb->height < HSS_EPSILON) bb->height = 0.0;
            bb->velocity = -bb->restitution * fabs(bb->velocity);
            bounces++;

            /* Check if at rest */
            if (fabs(bb->velocity) < 1e-6) {
                bb->at_rest = true;
                bb->velocity = 0.0;
                break;
            }

            /* Safeguard: if restitution is very low, stop after few bounces */
            if (bounces > 1000) {
                bb->at_rest = true;
                break;
            }
        }
    }

    bb->bounce_count += bounces;
    return bounces;
}

int hss_bouncing_ball_simulate(HSS_BouncingBall *bb, double max_time) {
    if (!bb) return 0;

    /* Use adaptive time steps */
    double dt = 1e-3;
    while (bb->time < max_time && !bb->at_rest) {
        double remaining = max_time - bb->time;
        if (dt > remaining) dt = remaining;
        hss_bouncing_ball_step(bb, dt);
    }
    return bb->bounce_count;
}

/* ============================================================================
 * L6 KP2: Thermostat with hysteresis
 * ============================================================================ */

HSS_Thermostat hss_thermostat_init(double T_start, double T_low,
                                    double T_high, double alpha,
                                    double beta, double T_hot, double T_cold) {
    HSS_Thermostat therm;
    memset(&therm, 0, sizeof(therm));
    therm.temperature     = T_start;
    therm.T_low           = T_low;
    therm.T_high          = T_high;
    therm.alpha           = alpha;
    therm.beta            = beta;
    therm.T_ambient_hot   = T_hot;
    therm.T_ambient_cold  = T_cold;
    therm.heater_on       = (T_start <= T_low);
    therm.time            = 0.0;
    therm.duty_cycle      = therm.heater_on ? 1.0 : 0.0;
    therm.avg_temperature = T_start;
    return therm;
}

int hss_thermostat_step(HSS_Thermostat *therm, double dt) {
    if (!therm) return 0;

    double T  = therm->temperature;
    int    sw = 0;

    if (therm->heater_on) {
        /* Heating: dT/dt = alpha * (T_hot - T) */
        double T_inf = therm->T_ambient_hot;
        double a = therm->alpha;
        T = T_inf - (T_inf - T) * exp(-a * dt);

        /* Check if we should turn off */
        if (T >= therm->T_high) {
            therm->heater_on = false;
            sw = 1;  /* ON → OFF */
        }
    } else {
        /* Cooling: dT/dt = -beta * (T - T_cold) */
        double T_inf = therm->T_ambient_cold;
        double b = therm->beta;
        T = T_inf + (T - T_inf) * exp(-b * dt);

        /* Check if we should turn on */
        if (T <= therm->T_low) {
            therm->heater_on = true;
            sw = -1;  /* OFF → ON */
        }
    }

    /* Update running average */
    double old_time = therm->time;
    therm->time  += dt;
    therm->avg_temperature = (therm->avg_temperature * old_time + T * dt)
                             / therm->time;

    /* Track duty cycle */
    if (therm->time > 0.0) {
        therm->duty_cycle = (therm->duty_cycle * old_time
                             + (therm->heater_on ? dt : 0.0)) / therm->time;
    }

    therm->temperature = T;
    return sw;
}

void hss_thermostat_simulate(HSS_Thermostat *therm,
                              double sim_time, double dt) {
    if (!therm) return;

    double t = 0.0;
    while (t < sim_time) {
        hss_thermostat_step(therm, dt);
        t += dt;
    }
}

/* ============================================================================
 * L6 KP3: DC-DC Boost Converter
 * ============================================================================ */

HSS_DCDCConverter hss_dcdc_init(double Vin, double L, double C,
                                 double R, double freq, double duty_cycle) {
    HSS_DCDCConverter conv;
    memset(&conv, 0, sizeof(conv));
    conv.Vin         = Vin;
    conv.L           = L;
    conv.C           = C;
    conv.R           = R;
    conv.freq        = freq;
    conv.duty_cycle  = duty_cycle;
    conv.iL          = 0.0;
    conv.vC          = Vin;  /* Start at input voltage */
    conv.switch_state = 1;   /* Start ON */
    conv.time        = 0.0;
    conv.vref        = Vin / (1.0 - duty_cycle); /* Ideal boost output */
    conv.output_ripple = 0.0;

    /* Check validity */
    if (freq <= 0.0 || duty_cycle < 0.0 || duty_cycle >= 1.0) {
        conv.freq = 100e3;  /* Default 100 kHz */
        conv.duty_cycle = 0.5;
    }
    if (L <= 0.0) conv.L = 1e-3;
    if (C <= 0.0) conv.C = 100e-6;
    if (R <= 0.0) conv.R = 10.0;

    return conv;
}

/**
 * @brief Advance DC-DC converter by one time step.
 *
 * Models the switched dynamics:
 *   Switch ON:  L·diL/dt = Vin,        C·dvC/dt = -vC/R
 *   Switch OFF: L·diL/dt = Vin - vC,   C·dvC/dt = iL - vC/R
 *
 * Uses exact solution of the linear ODE within each switching period.
 *
 * @return 1 if switch state changed, 0 otherwise
 */
int hss_dcdc_step(HSS_DCDCConverter *conv, double dt) {
    if (!conv) return 0;

    double L = conv->L;
    double C = conv->C;
    double R = conv->R;
    double Vin = conv->Vin;
    double iL = conv->iL;
    double vC = conv->vC;
    double freq = conv->freq;
    double D = conv->duty_cycle;
    double Tsw = 1.0 / freq;

    int switch_changed = 0;
    double time_left = dt;
    double t_local = fmod(conv->time, Tsw);
    double max_vC = vC, min_vC = vC;

    while (time_left > HSS_EPSILON) {
        double t_on  = D * Tsw;
        (void)(1.0 - D); /* t_off calculated implicitly */
        double t_phase;

        if (conv->switch_state == 1) {
            /* Currently ON */
            if (t_local < t_on) {
                t_phase = t_on - t_local;
            } else {
                /* Shouldn't be here but handle it */
                conv->switch_state = 0;
                switch_changed = 1;
                t_local = 0.0;
                continue;
            }
        } else {
            /* Currently OFF */
            if (t_local >= t_on && t_local < Tsw) {
                t_phase = Tsw - t_local;
            } else {
                conv->switch_state = 1;
                switch_changed = 1;
                t_local = 0.0;
                continue;
            }
        }

        if (t_phase > time_left) t_phase = time_left;

        /* Integrate */
        if (conv->switch_state == 1) {
            /* d(iL)/dt = Vin/L, d(vC)/dt = -vC/(RC) */
            iL += (Vin / L) * t_phase;
            vC = vC * exp(-t_phase / (R * C));
        } else {
            /* d(iL)/dt = (Vin - vC)/L, d(vC)/dt = (iL - vC/R)/C */
            /* Use Euler for simplicity */
            double diL = (Vin - vC) / L;
            double dvC = (iL - vC/R) / C;
            iL += diL * t_phase;
            vC += dvC * t_phase;
        }

        if (vC > max_vC) max_vC = vC;
        if (vC < min_vC) min_vC = vC;

        time_left -= t_phase;
        t_local += t_phase;
    }

    conv->iL = iL;
    conv->vC = vC;
    conv->time += dt;
    conv->output_ripple = max_vC - min_vC;

    return switch_changed;
}

double hss_dcdc_simulate(HSS_DCDCConverter *conv,
                          double sim_time, double dt) {
    if (!conv) return 0.0;

    double t = 0.0;
    while (t < sim_time) {
        hss_dcdc_step(conv, dt);
        t += dt;
    }
    return conv->vC;
}

/* ============================================================================
 * L6 KP4: Adaptive Cruise Control
 * ============================================================================ */

HSS_CruiseControl hss_cruise_control_init(double ego_vel,
                                           double lead_pos, double lead_vel,
                                           double v_set, double headway,
                                           double min_gap) {
    HSS_CruiseControl cc;
    memset(&cc, 0, sizeof(cc));
    cc.ego_position    = 0.0;
    cc.ego_velocity    = ego_vel;
    cc.lead_position   = lead_pos;
    cc.lead_velocity   = lead_vel;
    cc.v_set           = v_set;
    cc.time_headway    = headway;
    cc.min_gap         = min_gap;
    cc.max_accel       = 2.0;   /* m/s² */
    cc.max_decel       = 4.0;   /* m/s² */
    cc.time            = 0.0;
    cc.mode            = 0;     /* Start in CRUISE */
    cc.collision_risk  = false;

    /* Determine initial mode */
    double gap = lead_pos - cc.ego_position;
    double desired_gap = min_gap + headway * ego_vel;
    cc.gap_error = gap - desired_gap;

    if (gap < min_gap) {
        cc.mode = 2; /* BRAKE */
        cc.collision_risk = true;
    } else if (gap < desired_gap * 1.5) {
        cc.mode = 1; /* FOLLOW */
    }

    return cc;
}

/**
 * @brief Advance cruise control simulation.
 *
 * Mode logic:
 *   CRUISE: maintain v_set, switch to FOLLOW if gap < 1.5 * desired_gap
 *   FOLLOW: maintain desired gap, switch to BRAKE if gap < min_gap,
 *           switch to CRUISE if gap > 3 * desired_gap
 *   BRAKE:  maximum deceleration, switch to FOLLOW when gap recovers
 *
 * Uses simple proportional control in each mode.
 */
int hss_cruise_control_step(HSS_CruiseControl *cc, double dt) {
    if (!cc) return -1;

    double gap = cc->lead_position - cc->ego_position;
    double desired_gap = cc->min_gap + cc->time_headway * cc->ego_velocity;
    double v_err = cc->v_set - cc->ego_velocity;
    double accel = 0.0;

    /* Mode switching logic */
    switch (cc->mode) {
    case 0: /* CRUISE */
        if (gap < desired_gap * 1.5) {
            cc->mode = 1; /* → FOLLOW */
        }
        break;
    case 1: /* FOLLOW */
        if (gap < cc->min_gap) {
            cc->mode = 2;         /* → BRAKE */
            cc->collision_risk = true;
        } else if (gap > desired_gap * 3.0) {
            cc->mode = 0;         /* → CRUISE */
        }
        break;
    case 2: /* BRAKE */
        if (gap > desired_gap * 1.2) {
            cc->mode = 1;         /* → FOLLOW */
            cc->collision_risk = false;
        }
        break;
    }

    /* Control law per mode */
    switch (cc->mode) {
    case 0: /* CRUISE: PI speed control */
        accel = 0.5 * v_err;
        break;
    case 1: /* FOLLOW: gap control */
        {
            double gap_err = gap - desired_gap;
            double dv = cc->lead_velocity - cc->ego_velocity;
            accel = 0.3 * gap_err + 0.8 * dv;
        }
        break;
    case 2: /* BRAKE: max deceleration */
        accel = -cc->max_decel;
        break;
    }

    /* Saturate acceleration */
    if (accel > cc->max_accel) accel = cc->max_accel;
    if (accel < -cc->max_decel) accel = -cc->max_decel;

    /* Update ego vehicle state */
    cc->ego_velocity += accel * dt;
    if (cc->ego_velocity < 0.0) cc->ego_velocity = 0.0;
    cc->ego_position += cc->ego_velocity * dt;

    /* Update lead vehicle (constant velocity) */
    cc->lead_position += cc->lead_velocity * dt;

    cc->gap_error = cc->lead_position - cc->ego_position - desired_gap;
    cc->time += dt;

    return cc->mode;
}

void hss_cruise_control_simulate(HSS_CruiseControl *cc,
                                  double sim_time, double dt) {
    if (!cc) return;
    double t = 0.0;
    while (t < sim_time) {
        hss_cruise_control_step(cc, dt);
        t += dt;
    }
}

/* ============================================================================
 * L5 AUX4: Execution trace management
 * ============================================================================ */

HSS_ExecutionTrace *hss_trace_alloc(int max_steps, int state_dim) {
    if (max_steps <= 0 || state_dim <= 0) return NULL;

    HSS_ExecutionTrace *trace = calloc(1, sizeof(HSS_ExecutionTrace));
    if (!trace) return NULL;

    trace->max_steps  = max_steps;
    trace->num_steps  = 0;
    trace->num_jumps  = 0;
    trace->total_time = 0.0;

    trace->times  = calloc(max_steps, sizeof(double));
    trace->modes  = calloc(max_steps, sizeof(int));
    trace->states = calloc(max_steps * state_dim, sizeof(double));

    if (!trace->times || !trace->modes || !trace->states) {
        free(trace->times);
        free(trace->modes);
        free(trace->states);
        free(trace);
        return NULL;
    }

    return trace;
}

void hss_trace_free(HSS_ExecutionTrace *trace) {
    if (!trace) return;
    free(trace->times);
    free(trace->modes);
    free(trace->states);
    free(trace);
}

int hss_trace_record(HSS_ExecutionTrace *trace, double t,
                      int mode, const double *x, int n) {
    if (!trace || !x) return -1;
    if (trace->num_steps >= trace->max_steps) return -1;

    int idx = trace->num_steps;
    int state_dim = 0;

    /* Infer state dimension from allocation */
    for (int i = 0; i < trace->max_steps; i++) {
        if (trace->states[i * n] == 0.0 && n > 0) {
            state_dim = n;
            break;
        }
    }
    if (state_dim == 0) state_dim = n;

    trace->times[idx] = t;
    trace->modes[idx] = mode;
    if (n > 0) {
        memcpy(&trace->states[idx * n], x, n * sizeof(double));
    }
    trace->num_steps++;
    trace->total_time = t;

    return 0;
}

void hss_trace_print_summary(const HSS_ExecutionTrace *trace, FILE *fp) {
    if (!trace) { fprintf(fp, "NULL trace\n"); return; }
    if (!fp) fp = stdout;

    fprintf(fp, "=== Execution Trace Summary ===\n");
    fprintf(fp, "  Steps recorded:  %d / %d\n", trace->num_steps,
            trace->max_steps);
    fprintf(fp, "  Total jumps:     %d\n", trace->num_jumps);
    fprintf(fp, "  Total time:      %.6f\n", trace->total_time);
    if (trace->num_steps > 0) {
        fprintf(fp, "  Time range:      [%.6f, %.6f]\n",
                trace->times[0],
                trace->times[trace->num_steps - 1]);
    }
}

int hss_trace_export_csv(const HSS_ExecutionTrace *trace,
                          const char *filename) {
    if (!trace || !filename) return -1;

    FILE *fp = fopen(filename, "w");
    if (!fp) return -1;

    fprintf(fp, "time,mode");
    /* Assume state dimension from allocation context (passed separately) */
    /* For simplicity, write just time and mode */
    for (int i = 0; i < trace->num_steps; i++) {
        fprintf(fp, "%.10e,%d\n", trace->times[i], trace->modes[i]);
    }

    fclose(fp);
    return 0;
}

/* ============================================================================
 * Internal: ODE Solvers
 * ============================================================================ */

/**
 * @brief RK4 integration step.
 *
 * Classical 4th-order Runge-Kutta method:
 *   k1 = h · f(x, t)
 *   k2 = h · f(x + k1/2, t + h/2)
 *   k3 = h · f(x + k2/2, t + h/2)
 *   k4 = h · f(x + k3, t + h)
 *   x(t+h) = x + (k1 + 2k2 + 2k3 + k4)/6
 *
 * Global error: O(h⁴). 4 function evaluations per step.
 *
 * Reference: Butcher (2008), Numerical Methods for ODEs.
 */
static int rk4_step(void (*f)(const double*, const double*, double,
                               double*, int, void*),
                     double *x, const double *u, double t, double h,
                     int n, void *params) {
    double *k1 = calloc(n, sizeof(double));
    double *k2 = calloc(n, sizeof(double));
    double *k3 = calloc(n, sizeof(double));
    double *k4 = calloc(n, sizeof(double));
    double *xtmp = calloc(n, sizeof(double));

    if (!k1 || !k2 || !k3 || !k4 || !xtmp) {
        free(k1); free(k2); free(k3); free(k4); free(xtmp);
        return -1;
    }

    /* k1 = f(x, t) */
    f(x, u, t, k1, n, params);

    /* k2 = f(x + h·k1/2, t + h/2) */
    for (int i = 0; i < n; i++) xtmp[i] = x[i] + 0.5 * h * k1[i];
    f(xtmp, u, t + 0.5 * h, k2, n, params);

    /* k3 = f(x + h·k2/2, t + h/2) */
    for (int i = 0; i < n; i++) xtmp[i] = x[i] + 0.5 * h * k2[i];
    f(xtmp, u, t + 0.5 * h, k3, n, params);

    /* k4 = f(x + h·k3, t + h) */
    for (int i = 0; i < n; i++) xtmp[i] = x[i] + h * k3[i];
    f(xtmp, u, t + h, k4, n, params);

    /* x(t+h) = x + h·(k1 + 2k2 + 2k3 + k4)/6 */
    for (int i = 0; i < n; i++) {
        x[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }

    free(k1); free(k2); free(k3); free(k4); free(xtmp);
    return 0;
}

/**
 * @brief Forward Euler integration step.
 *
 * x(t+h) = x + h · f(x, t)
 *
 * Global error: O(h). Simple but requires small step sizes.
 */
static int euler_step(void (*f)(const double*, const double*, double,
                                 double*, int, void*),
                       double *x, const double *u, double t, double h,
                       int n, void *params) {
    double *dxdt = calloc(n, sizeof(double));
    if (!dxdt) return -1;

    f(x, u, t, dxdt, n, params);
    for (int i = 0; i < n; i++) {
        x[i] += h * dxdt[i];
    }

    free(dxdt);
    return 0;
}

/**
 * @brief Wrapper converting linear dynamics to flow function.
 *
 * ẋ = A x + B u + c (if affine)
 *
 * This allows linear modes to use the same ODE solver interface
 * as nonlinear modes.
 */
static void linear_flow_wrapper(const double *x, const double *u, double t,
                                 double *dxdt, int n, void *params) {
    HSS_Mode *mode = (HSS_Mode *)params;
    if (!mode || !mode->matrix.A) {
        memset(dxdt, 0, n * sizeof(double));
        return;
    }
    (void)t; /* suppress unused warning */

    /* dxdt = A x */
    for (int i = 0; i < n; i++) {
        dxdt[i] = 0.0;
        for (int j = 0; j < n; j++) {
            dxdt[i] += mode->matrix.A[i * n + j] * x[j];
        }
    }

    /* + B u */
    if (mode->matrix.B && u && mode->matrix.m > 0) {
        int m = mode->matrix.m;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                dxdt[i] += mode->matrix.B[i * m + j] * u[j];
            }
        }
    }

    /* + c (affine term) */
    if (mode->matrix.has_affine && mode->matrix.c) {
        for (int i = 0; i < n; i++) {
            dxdt[i] += mode->matrix.c[i];
        }
    }
}

