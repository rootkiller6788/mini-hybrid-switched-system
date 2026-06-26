#ifndef SWITCHED_APPLICATIONS_H
#define SWITCHED_APPLICATIONS_H

#include "switched_types.h"

/* ============================================================================
 * Switched System Applications (L6, L7)
 *
 * Domain-specific implementations of switched system models:
 *   - DC-DC power converters (power electronics)
 *   - Thermostat control (HVAC/building automation)
 *   - Vehicle spacing control (automotive/transportation)
 *   - Networked control with dropouts (NCS/IoT)
 * ============================================================================ */

/* --------------------------------------------------------------------------
 * L6: Canonical Problem - DC-DC Boost Converter
 * -------------------------------------------------------------------------- */

/**
 * Initialize a DC-DC boost converter model.
 *
 * Boost converter: steps up input voltage Vin to higher output Vout.
 * State: x = [iL, vC]^T (inductor current, capacitor voltage).
 *
 * @param Vin     Input voltage (e.g., 12V)
 * @param Vout    Target output voltage
 * @param L       Inductance in Henry
 * @param C       Capacitance in Farad
 * @param R_load  Load resistance in Ohm
 * @param freq    Switching frequency in Hz
 * @return        Initialized DCDC converter
 */
DCDCConverter* dcdc_create(double Vin, double Vout, double L, double C,
                            double R_load, double freq);
void dcdc_free(DCDCConverter *conv);

/**
 * Compute the A matrices for ON and OFF modes of the boost converter.
 *
 * ON mode (switch closed):  diL/dt = Vin/L,            dvC/dt = -vC/(R*C)
 * OFF mode (switch open):   diL/dt = (Vin - vC)/L,     dvC/dt = (iL - vC/R)/C
 *
 * @param conv   Converter parameters
 * @param A_on   Output: A matrix for ON mode (2x2)
 * @param A_off  Output: A matrix for OFF mode (2x2)
 */
void dcdc_get_matrices(const DCDCConverter *conv, SwitchedMatrix *A_on, SwitchedMatrix *A_off);

/**
 * Compute the equilibrium (average) state for a given duty cycle D.
 *
 * For a boost converter in continuous conduction mode:
 *   Vout = Vin / (1 - D)
 *   iL_ss = Vin / (R * (1-D)^2)
 *
 * @param conv  Converter parameters
 * @param D     Duty cycle D in (0, 1)
 * @param iL_ss Output: steady-state inductor current
 * @param vC_ss Output: steady-state capacitor voltage
 */
void dcdc_equilibrium(const DCDCConverter *conv, double D, double *iL_ss, double *vC_ss);

/**
 * Simulate DC-DC converter with PWM switching.
 *
 * PWM pattern: ON for D*T, OFF for (1-D)*T where T = 1/freq.
 * Uses exact discretization for each mode.
 *
 * @param conv   Converter model
 * @param t_end  Simulation duration
 * @param dt     Time step (should be << T for accuracy)
 */
void dcdc_simulate(DCDCConverter *conv, double t_end, double dt);

/**
 * Analyze stability of the DC-DC converter as a switched system.
 * Computes stability margin and verifies dwell-time conditions.
 */
bool dcdc_analyze_stability(const DCDCConverter *conv);

/**
 * Compute ripple (peak-to-peak variation) in steady state.
 * Ripple in vC is approximately: delta_vC = (Vout * D) / (R * C * freq)
 */
void dcdc_compute_ripple(const DCDCConverter *conv, double D, double *ripple_iL, double *ripple_vC);

/* --------------------------------------------------------------------------
 * L6: Canonical Problem - Thermostat Control
 * -------------------------------------------------------------------------- */

/**
 * Initialize a thermostat (bang-bang) control system.
 *
 * Hysteresis-based switching between OFF, HEATING, COOLING modes.
 * Switching logic:
 *   - temp < setpoint - deadband: switch to HEATING
 *   - temp > setpoint + deadband: switch to COOLING
 *   - otherwise: OFF
 *
 * @param setpoint      Desired temperature
 * @param deadband      Hysteresis width (+/- deadband)
 * @param heating_rate  Temperature increase per second (heating)
 * @param cooling_rate  Temperature decrease per second (cooling)
 * @param ambient_loss  Natural heat loss rate
 * @return              Initialized thermostat system
 */
ThermostatSystem* thermo_create(double setpoint, double deadband,
                                 double heating_rate, double cooling_rate,
                                 double ambient_loss);
void thermo_free(ThermostatSystem *thermo);

/**
 * Update thermostat state by one time step.
 * Determines the active mode based on current temperature and
 * advances the dynamics accordingly.
 *
 * @param thermo  Thermostat system
 * @param dt      Time step
 */
void thermo_step(ThermostatSystem *thermo, double dt);

/**
 * Run thermostat simulation and compute switching statistics.
 * Reports number of mode switches, average ON time, duty cycle.
 */
void thermo_simulate(ThermostatSystem *thermo, double t_end, double dt);

/**
 * Analyze thermostat as a switched system: compute equivalent
 * A matrices for each mode and assess stability.
 */
void thermo_switched_analysis(const ThermostatSystem *thermo);

/* --------------------------------------------------------------------------
 * L6: Vehicle Spacing Control
 * -------------------------------------------------------------------------- */

/**
 * Initialize vehicle spacing control with 4 operating modes:
 *   0=CRUISE (maintain speed), 1=ACCEL (increase speed),
 *   2=DECEL (reduce speed), 3=EMERGENCY (hard brake).
 *
 * @param ego_speed   Initial ego vehicle speed (m/s)
 * @param lead_speed  Lead vehicle speed (m/s)
 * @param gap         Initial gap (m)
 * @param safe_gap    Desired safe gap (m)
 * @return            Initialized vehicle spacing controller
 */
VehicleSpacingControl* vsc_create(double ego_speed, double lead_speed,
                                   double gap, double safe_gap);
void vsc_free(VehicleSpacingControl *vsc);

/**
 * Determine the active mode based on gap error and relative speed.
 * Mode selection uses hysteresis to prevent chattering.
 *
 * @param vsc  Vehicle spacing controller
 * @return     Active mode index
 */
int vsc_determine_mode(const VehicleSpacingControl *vsc);

/**
 * Update vehicle state based on current mode.
 * Each mode has different acceleration limits and control laws.
 */
void vsc_step(VehicleSpacingControl *vsc, double dt);

/**
 * Simulate vehicle spacing over time and assess stability
 * as a switched system under different lead vehicle profiles.
 */
void vsc_simulate(VehicleSpacingControl *vsc, double t_end, double dt);

/* --------------------------------------------------------------------------
 * L7: Networked Control with Packet Dropouts
 * -------------------------------------------------------------------------- */

/**
 * Initialize networked control system with dropout model.
 *
 * @param state_dim         Plant state dimension
 * @param ctrl_gain         Controller gain matrix (flattened, state_dim^2)
 * @param dropout_rate      Packet loss probability [0, 1]
 * @param max_allowable_loss  MADB (Maximum Allowable Dropout Bound)
 * @return                  Initialized NCS dropout model
 */
NetworkedControlDropout* ncs_create(int state_dim, const double *ctrl_gain,
                                     double dropout_rate, double max_allowable_loss);
void ncs_free(NetworkedControlDropout *ncs);

/**
 * Simulate one step of networked control with random dropout.
 * If packet is lost, the previous control value is held (zero-order hold).
 *
 * @param ncs  NCS dropout model
 * @param dt   Time step
 */
void ncs_step(NetworkedControlDropout *ncs, double dt);

/**
 * Compute the MADB (Maximum Allowable Dropout Bound).
 *
 * For a system x_{k+1} = A x_k + B u_k with u_k = K x_k:
 * If n consecutive packets are lost, u_k = u_{k-n} (hold strategy).
 *
 * MADB is the maximum n for which the switched system
 * (between closed-loop and open-loop modes) remains stable.
 *
 * @param A      Plant A matrix
 * @param B      Plant B matrix
 * @param K      Controller gain
 * @param n      State dimension
 * @return       Maximum allowable consecutive dropouts
 */
int ncs_compute_madb(const SwitchedMatrix *A, const SwitchedMatrix *B,
                     const SwitchedMatrix *K, int n);

/**
 * Simulate NCS over extended time and report statistics:
 * - Total packet loss count
 * - Maximum consecutive losses
 * - Stability flag
 */
void ncs_simulate(NetworkedControlDropout *ncs, double t_end, double dt);

/**
 * Analyze the switched system stability of the NCS under dropout.
 * The system switches between:
 *   - Mode 0: Packet received (closed-loop, stable)
 *   - Mode 1: Packet lost (open-loop, may be unstable)
 *
 * Uses average dwell-time analysis to determine if the dropout
 * pattern allows GUES.
 */
bool ncs_stability_analysis(NetworkedControlDropout *ncs);

#endif /* SWITCHED_APPLICATIONS_H */
