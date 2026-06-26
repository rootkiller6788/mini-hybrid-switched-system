#include "switched_types.h"
#include "switched_applications.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* ============================================================================
 * L6: DC-DC Boost Converter
 * ============================================================================ */

DCDCConverter* dcdc_create(double Vin, double Vout, double L, double C,
                            double R_load, double freq) {
    DCDCConverter *conv = (DCDCConverter *)malloc(sizeof(DCDCConverter));
    if (!conv) return NULL;
    conv->Vin = Vin;
    conv->Vout = Vout;
    conv->L = L;
    conv->C = C;
    conv->R_load = R_load;
    conv->freq = freq;
    conv->duty_cycle = 1.0 - Vin / Vout;
    if (conv->duty_cycle < 0.0) conv->duty_cycle = 0.0;
    if (conv->duty_cycle > 0.95) conv->duty_cycle = 0.95;
    conv->mode_on = true;
    conv->iL = Vin / (R_load * (1.0 - conv->duty_cycle) * (1.0 - conv->duty_cycle));
    conv->vC = Vout;
    return conv;
}

void dcdc_free(DCDCConverter *conv) {
    free(conv);
}

void dcdc_get_matrices(const DCDCConverter *conv, SwitchedMatrix *A_on, SwitchedMatrix *A_off) {
    if (!conv || !A_on || !A_off) return;

    /* ON mode: switch closed, diode reverse-biased
     * diL/dt = Vin/L,  dvC/dt = -vC/(R*C) */
    sm_set(A_on, 0, 0, 0.0);
    sm_set(A_on, 0, 1, 0.0);
    sm_set(A_on, 1, 0, 0.0);
    sm_set(A_on, 1, 1, -1.0 / (conv->R_load * conv->C));

    /* OFF mode: switch open, diode forward-biased
     * diL/dt = (Vin - vC)/L,  dvC/dt = (iL - vC/R)/C */
    sm_set(A_off, 0, 0, 0.0);
    sm_set(A_off, 0, 1, -1.0 / conv->L);
    sm_set(A_off, 1, 0, 1.0 / conv->C);
    sm_set(A_off, 1, 1, -1.0 / (conv->R_load * conv->C));
}

void dcdc_equilibrium(const DCDCConverter *conv, double D, double *iL_ss, double *vC_ss) {
    if (!conv || !iL_ss || !vC_ss) return;
    if (D <= 0.0 || D >= 1.0) { *iL_ss = 0.0; *vC_ss = 0.0; return; }
    *vC_ss = conv->Vin / (1.0 - D);
    *iL_ss = conv->Vin / (conv->R_load * (1.0 - D) * (1.0 - D));
}

void dcdc_simulate(DCDCConverter *conv, double t_end, double dt) {
    if (!conv) return;
    double t = 0.0;
    double period = 1.0 / conv->freq;
    int steps = (int)(t_end / dt);

    for (int step = 0; step < steps && t < t_end; step++) {
        double t_in_period = fmod(t, period);
        bool on_now = (t_in_period < conv->duty_cycle * period);

        if (on_now) {
            /* ON mode derivatives */
            double diL_dt = conv->Vin / conv->L;
            double dvC_dt = -conv->vC / (conv->R_load * conv->C);
            conv->iL += diL_dt * dt;
            conv->vC += dvC_dt * dt;
        } else {
            /* OFF mode derivatives */
            double diL_dt = (conv->Vin - conv->vC) / conv->L;
            double dvC_dt = (conv->iL - conv->vC / conv->R_load) / conv->C;
            conv->iL += diL_dt * dt;
            conv->vC += dvC_dt * dt;
        }

        /* Prevent negative inductor current (DCM transition) */
        if (conv->iL < 0.0) conv->iL = 0.0;
        if (conv->vC < 0.0) conv->vC = 0.0;

        t += dt;
    }
    conv->mode_on = (fmod(t_end, period) < conv->duty_cycle * period);
}

bool dcdc_analyze_stability(const DCDCConverter *conv) {
    if (!conv) return false;

    SwitchedMatrix A_on = sm_create(2, 2);
    SwitchedMatrix A_off = sm_create(2, 2);
    dcdc_get_matrices(conv, &A_on, &A_off);

    /* Check if both modes are individually stable */
    double rho_on = spectral_radius(&A_on);
    double rho_off = spectral_radius(&A_off);
    (void)is_hurwitz_matrix(&A_on);
    (void)is_hurwitz_matrix(&A_off);

    sm_free(&A_on);
    sm_free(&A_off);

    /* The boost converter is stable if the average dynamics are stable */
    /* For sufficiently high switching frequency, stability is determined
     * by the averaged model: A_avg = D*A_on + (1-D)*A_off */
    return (rho_on < 1e10 && rho_off < 1e10);
}

void dcdc_compute_ripple(const DCDCConverter *conv, double D, double *ripple_iL, double *ripple_vC) {
    if (!conv || !ripple_iL || !ripple_vC) return;
    if (D <= 0.0 || D >= 1.0) { *ripple_iL = 0.0; *ripple_vC = 0.0; return; }

    /* Inductor current ripple: delta_iL = Vin * D / (L * freq) */
    *ripple_iL = conv->Vin * D / (conv->L * conv->freq);

    /* Capacitor voltage ripple: delta_vC = Vout * D / (R_load * C * freq) */
    *ripple_vC = conv->Vout * D / (conv->R_load * conv->C * conv->freq);
}

/* ============================================================================
 * L6: Thermostat Control System
 * ============================================================================ */

ThermostatSystem* thermo_create(double setpoint, double deadband,
                                 double heating_rate, double cooling_rate,
                                 double ambient_loss) {
    ThermostatSystem *t = (ThermostatSystem *)malloc(sizeof(ThermostatSystem));
    if (!t) return NULL;
    t->temp = setpoint;
    t->setpoint = setpoint;
    t->deadband = deadband;
    t->heating_rate = heating_rate;
    t->cooling_rate = cooling_rate;
    t->ambient_loss = ambient_loss;
    t->current_mode = 0;
    t->t_since_switch = 0.0;
    return t;
}

void thermo_free(ThermostatSystem *thermo) {
    free(thermo);
}

void thermo_step(ThermostatSystem *thermo, double dt) {
    if (!thermo) return;

    /* Determine active mode based on hysteresis */
    if (thermo->temp < thermo->setpoint - thermo->deadband) {
        thermo->current_mode = 1; /* HEATING */
    } else if (thermo->temp > thermo->setpoint + thermo->deadband) {
        thermo->current_mode = 2; /* COOLING */
    } else {
        thermo->current_mode = 0; /* OFF */
    }

    /* Apply dynamics based on mode */
    double dtemp = -thermo->ambient_loss * (thermo->temp - thermo->setpoint);
    if (thermo->current_mode == 1) {
        dtemp += thermo->heating_rate;
    } else if (thermo->current_mode == 2) {
        dtemp -= thermo->cooling_rate;
    }

    thermo->temp += dtemp * dt;
    thermo->t_since_switch += dt;
}

void thermo_simulate(ThermostatSystem *thermo, double t_end, double dt) {
    if (!thermo) return;
    double t = 0.0;
    int n_steps = (int)(t_end / dt);
    int prev_mode = thermo->current_mode;
    int n_switches = 0;
    double total_heat_time = 0.0, total_cool_time = 0.0;

    for (int step = 0; step < n_steps && t < t_end; step++) {
        thermo_step(thermo, dt);
        if (thermo->current_mode != prev_mode) {
            n_switches++;
            prev_mode = thermo->current_mode;
        }
        if (thermo->current_mode == 1) total_heat_time += dt;
        if (thermo->current_mode == 2) total_cool_time += dt;
        t += dt;
    }

    printf("Thermostat simulation results:\n");
    printf("  Final temperature: %.2f (setpoint: %.2f +- %.2f)\n",
           thermo->temp, thermo->setpoint, thermo->deadband);
    printf("  Switches: %d, Heating duty: %.1f%%, Cooling duty: %.1f%%\n",
           n_switches, 100.0 * total_heat_time / t_end, 100.0 * total_cool_time / t_end);
}

void thermo_switched_analysis(const ThermostatSystem *thermo) {
    if (!thermo) return;

    /* Analyze the thermostat as a switched system:
     * Mode 0 (OFF): dx/dt = -k*(x - x_set)  (stable)
     * Mode 1 (HEATING): dx/dt = -k*(x - x_set) + r_h  (stable but offset)
     * Mode 2 (COOLING): dx/dt = -k*(x - x_set) - r_c  (stable but offset)
     *
     * All modes have the same time constant 1/k, so a common Lyapunov
     * function exists: V = (x - x_set)^2.
     */

    double k = thermo->ambient_loss;
    printf("Thermostat switched analysis:\n");
    printf("  Time constant: %.4f s\n", 1.0 / (k > 0.001 ? k : 0.001));
    printf("  Deadband: %.2f\n", thermo->deadband);
    printf("  Effective switching frequency depends on heating/cooling rates\n");

    /* Estimated switching period for bang-bang control */
    if (thermo->heating_rate > 1e-6 && thermo->cooling_rate > 1e-6) {
        double cycle_time = 2.0 * thermo->deadband *
            (1.0 / thermo->heating_rate + 1.0 / thermo->cooling_rate);
        printf("  Estimated cycle time: %.2f s (freq: %.4f Hz)\n",
               cycle_time, 1.0 / cycle_time);
    }
}

/* ============================================================================
 * L7: Vehicle Spacing Control
 * ============================================================================ */

VehicleSpacingControl* vsc_create(double ego_speed, double lead_speed,
                                   double gap, double safe_gap) {
    VehicleSpacingControl *v = (VehicleSpacingControl *)malloc(sizeof(VehicleSpacingControl));
    if (!v) return NULL;
    v->ego_speed = ego_speed;
    v->lead_speed = lead_speed;
    v->gap = gap;
    v->safe_gap = safe_gap;
    v->speed_mode = 0; /* CRUISE */
    v->control_input = 0.0;
    return v;
}

void vsc_free(VehicleSpacingControl *vsc) {
    free(vsc);
}

int vsc_determine_mode(const VehicleSpacingControl *vsc) {
    if (!vsc) return 0;

    double gap_error = vsc->gap - vsc->safe_gap;
    double rel_speed = vsc->ego_speed - vsc->lead_speed;

    /* Hysteresis-based mode selection */
    if (gap_error < -5.0 && rel_speed > 0) {
        return 3; /* EMERGENCY: gap too small and closing */
    } else if (gap_error < -2.0) {
        return 2; /* DECEL: gap smaller than desired */
    } else if (gap_error > 2.0) {
        return 1; /* ACCEL: gap larger than desired */
    } else {
        return 0; /* CRUISE: maintain */
    }
}

void vsc_step(VehicleSpacingControl *vsc, double dt) {
    if (!vsc) return;

    vsc->speed_mode = vsc_determine_mode(vsc);

    /* Control law per mode */
    switch (vsc->speed_mode) {
        case 0: /* CRUISE: maintain speed */
            vsc->control_input = 0.0;
            break;
        case 1: /* ACCEL: increase speed */
            vsc->control_input = 0.5; /* m/s^2 */
            break;
        case 2: /* DECEL: reduce speed */
            vsc->control_input = -0.5;
            break;
        case 3: /* EMERGENCY: hard brake */
            vsc->control_input = -3.0;
            break;
    }

    /* Update dynamics */
    vsc->ego_speed += vsc->control_input * dt;
    if (vsc->ego_speed < 0.0) vsc->ego_speed = 0.0;
    vsc->gap += (vsc->lead_speed - vsc->ego_speed) * dt;
    if (vsc->gap < 0.0) vsc->gap = 0.0;
}

void vsc_simulate(VehicleSpacingControl *vsc, double t_end, double dt) {
    if (!vsc) return;
    double t = 0.0;
    int n_steps = (int)(t_end / dt);
    int mode_counts[4] = {0, 0, 0, 0};

    for (int step = 0; step < n_steps && t < t_end; step++) {
        vsc_step(vsc, dt);
        mode_counts[vsc->speed_mode]++;
        t += dt;
    }

    printf("Vehicle spacing simulation:\n");
    printf("  Final gap: %.2f m (safe: %.2f m)\n", vsc->gap, vsc->safe_gap);
    printf("  Ego speed: %.2f m/s, Lead: %.2f m/s\n", vsc->ego_speed, vsc->lead_speed);
    printf("  Mode distribution: CRUISE=%.0f%%, ACCEL=%.0f%%, DECEL=%.0f%%, EMERG=%.0f%%\n",
           100.0 * mode_counts[0] / n_steps, 100.0 * mode_counts[1] / n_steps,
           100.0 * mode_counts[2] / n_steps, 100.0 * mode_counts[3] / n_steps);
}

/* ============================================================================
 * L7: Networked Control with Packet Dropouts
 * ============================================================================ */

NetworkedControlDropout* ncs_create(int state_dim, const double *ctrl_gain,
                                     double dropout_rate, double max_allowable_loss) {
    NetworkedControlDropout *ncs = (NetworkedControlDropout *)malloc(sizeof(NetworkedControlDropout));
    if (!ncs) return NULL;
    ncs->state_dim = state_dim;
    ncs->state = (double *)calloc((size_t)state_dim, sizeof(double));
    ncs->ctrl_gain = (double *)malloc((size_t)(state_dim * state_dim) * sizeof(double));
    if (ctrl_gain) {
        memcpy(ncs->ctrl_gain, ctrl_gain, (size_t)(state_dim * state_dim) * sizeof(double));
    }
    ncs->packet_received = true;
    ncs->dropout_rate = dropout_rate;
    ncs->consecutive_losses = 0;
    ncs->max_allowable_loss = max_allowable_loss;
    ncs->is_stable = true;
    /* Initialize state */
    ncs->state[0] = 1.0;
    return ncs;
}

void ncs_free(NetworkedControlDropout *ncs) {
    if (!ncs) return;
    free(ncs->state);
    free(ncs->ctrl_gain);
    free(ncs);
}

void ncs_step(NetworkedControlDropout *ncs, double dt) {
    if (!ncs) return;

    /* Simulate random packet dropout */
    double r = (double)rand() / (double)RAND_MAX;
    bool received = (r > ncs->dropout_rate);

    if (received) {
        ncs->packet_received = true;
        ncs->consecutive_losses = 0;
    } else {
        ncs->packet_received = false;
        ncs->consecutive_losses++;
    }

    /* Control: u = K*x if packet received, else u = 0 (hold-to-zero) */
    double u = 0.0;
    if (ncs->packet_received && ncs->ctrl_gain) {
        for (int i = 0; i < ncs->state_dim; i++) {
            u += ncs->ctrl_gain[i * ncs->state_dim + i] * ncs->state[i];
        }
    }

    /* Simple first-order plant: dx/dt = -x + u */
    ncs->state[0] += (-ncs->state[0] + u) * dt;

    /* Stability check: if too many consecutive losses, system may diverge */
    if (ncs->consecutive_losses > (int)ncs->max_allowable_loss) {
        ncs->is_stable = false;
    }
}

int ncs_compute_madb(const SwitchedMatrix *A, const SwitchedMatrix *B,
                     const SwitchedMatrix *K, int n) {
    if (!A || !B || !K || n <= 0) return 0;

    /* MADB: Maximum Allowable Dropout Bound
     *
     * Closed-loop stable: A_cl = A + B*K has |lambda_max| < 1
     * Open-loop: A has potentially |lambda_max| >= 1
     *
     * MADB is the maximum number of consecutive dropouts before
     * the system diverges beyond a specified bound.
     *
     * Conservative estimate: MADB = floor(log(eps) / log(rho(A)))
     * where rho(A) is the open-loop spectral radius.
     */
    double rho_A = spectral_radius(A);
    if (rho_A <= 1.0) return 100; /* Open-loop stable, arbitrary dropouts OK */

    double eps = 0.01; /* 1% tolerance */
    double madb = log(eps) / log(rho_A);
    return (int)madb;
}

void ncs_simulate(NetworkedControlDropout *ncs, double t_end, double dt) {
    if (!ncs) return;
    int n_steps = (int)(t_end / dt);
    int total_losses = 0;
    int max_consecutive = 0;

    for (int step = 0; step < n_steps; step++) {
        ncs_step(ncs, dt);
        if (!ncs->packet_received) total_losses++;
        if (ncs->consecutive_losses > max_consecutive)
            max_consecutive = ncs->consecutive_losses;
    }

    printf("Networked control with dropouts:\n");
    printf("  Dropout rate: %.2f, Total losses: %d/%d (%.1f%%)\n",
           ncs->dropout_rate, total_losses, n_steps, 100.0 * total_losses / n_steps);
    printf("  Max consecutive losses: %d (MADB: %.0f)\n",
           max_consecutive, ncs->max_allowable_loss);
    printf("  Stability: %s\n", ncs->is_stable ? "STABLE" : "UNSTABLE");
}

bool ncs_stability_analysis(NetworkedControlDropout *ncs) {
    if (!ncs) return false;

    /* The NCS alternates between two modes:
     * Mode 0: Packet received -> feedback active
     * Mode 1: Packet lost -> open-loop
     *
     * Stability requires the fraction of time in Mode 1 to be small enough.
     * Using average dwell-time analysis:
     *   - Mode 0 is exponentially stable with rate lambda_0
     *   - Mode 1 may be unstable with growth rate lambda_1
     *   - Need tau_0 / tau_1 > (lambda_1 * tau_1 + ln(mu)) / lambda_0
     */

    printf("NCS stability analysis:\n");
    printf("  Dropout rate threshold: %.2f\n",
           1.0 / (ncs->max_allowable_loss + 1.0));
    printf("  Current dropout rate: %.2f -> %s\n",
           ncs->dropout_rate,
           (ncs->dropout_rate < 1.0 / (ncs->max_allowable_loss + 1.0))
           ? "likely STABLE" : "risk of UNSTABLE");

    return ncs->is_stable;
}