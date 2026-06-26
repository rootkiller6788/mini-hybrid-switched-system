/* reset_element.h - Reset Control Elements: Clegg, FORE, Reset PID
 *
 * Defines the canonical reset control elements: Clegg integrator,
 * First-Order Reset Element (FORE), Second-Order Reset Element (SORE),
 * Reset PID controller, and Reset Lead-Lag compensator.
 *
 * Knowledge coverage:
 *   L1: Clegg integrator, FORE, SORE, reset PID (definitions)
 *   L2: Describing function, sinusoidal input analysis, phase lead
 *       from reset, reset-induced fundamental component
 *   L3: Transfer function base + reset map composition
 *   L5: Step/impulse/sine response evaluation algorithms
 *
 * Ref: [BB12] Ch.3 for Clegg/FORE/SORE; [HR75] for FORE;
 *      [ZCH00] for reset PID; Hazeleger, Heertjes, Nijmeijer (2011)
 *      for second-order reset elements.
 */

#ifndef RESET_ELEMENT_H
#define RESET_ELEMENT_H

#include "reset_core.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L1: Reset Element Types
 * ================================================================ */

/** ResetElementType - classifies the canonical reset controller elements. */
typedef enum {
    RESET_ELEM_CLEGG        = 0,  /* Clegg integrator (1st order)      */
    RESET_ELEM_FORE         = 1,  /* First Order Reset Element         */
    RESET_ELEM_SORE         = 2,  /* Second Order Reset Element        */
    RESET_ELEM_RESET_PID    = 3,  /* Reset PID: P+I_reset+D            */
    RESET_ELEM_RESET_LEAD   = 4,  /* Reset lead compensator            */
    RESET_ELEM_RESET_LAG    = 5,  /* Reset lag compensator             */
    RESET_ELEM_CUSTOM       = 6   /* User-defined reset system         */
} ResetElementType;

/* ================================================================
 * L1: Clegg Integrator (Clegg, 1958)
 * ================================================================ */

/** CleggIntegrator - the foundational reset element.
 *
 *  Dynamics (between resets):
 *    dx/dt = e(t)              (pure integration)
 *    y(t)  = x(t)
 *
 *  Reset law (at zero crossing of e):
 *    x(t_k^+) = 0              (full state erasure)
 *
 *  Fundamental property: the describing function of the Clegg
 *  integrator provides 51.85 degrees of phase lag, compared to
 *  90 degrees for a linear integrator. This 38.15-degree phase
 *  advantage is the key benefit of reset control.
 *
 *  Describing Function:
 *    N_Clegg(jw) = (1.62 / w) * exp(-j * 51.85 deg)
 *    vs. N_Linear(jw) = (1.0 / w) * exp(-j * 90 deg)
 *
 *  Ref: [Cle58]; [BB12] Section 3.2.1
 */
typedef struct {
    ResetSystem    *base;       /**< underlying reset system           */
    double          x0;         /**< initial integrator state          */
    double          y;          /**< current output                    */
    int             n_resets;   /**< reset count for this element      */
} CleggIntegrator;

/* ================================================================
 * L1: FORE - First Order Reset Element (Horowitz & Rosenbaum, 1975)
 * ================================================================ */

/** ForeElement - First Order Reset Element.
 *
 *  Base transfer function:
 *    G(s) = K / (tau * s + 1)   (1st order low-pass filter)
 *
 *  State-space (flow):
 *    dx/dt = (-1/tau) * x + (K/tau) * e
 *    y     = x
 *
 *  Reset law:
 *    x(t_k^+) = rho * x(t_k^-),   rho in [0, 1)
 *
 *  The FORE generalizes the Clegg integrator: Clegg = FORE with
 *  tau -> infinity and rho = 0.
 *
 *  Describing Function:
 *    N_FORE(jw) = K * (1 + j*w*tau*rho) / ((1 + j*w*tau) * ...)
 *    Phase advantage depends on rho and frequency.
 *
 *  Ref: [HR75]; [BB12] Section 3.3
 */
typedef struct {
    ResetSystem    *base;       /**< underlying reset system           */
    double          K;          /**< DC gain                          */
    double          tau;        /**< time constant (seconds)          */
    double          rho;        /**< reset ratio                      */
    double          x0;         /**< initial state                    */
    double          y;          /**< current output                   */
} ForeElement;

/* ================================================================
 * L1: SORE - Second Order Reset Element
 * ================================================================ */

/** SoreElement - Second Order Reset Element.
 *
 *  Base transfer function:
 *    G(s) = wn^2 / (s^2 + 2*zeta*wn*s + wn^2)
 *
 *  State-space (controllable canonical form):
 *    A = [   0         1    ]
 *        [ -wn^2   -2*zeta*wn ]
 *    B = [ 0; wn^2 ]
 *    C = [ 1,  0 ]
 *    D = [ 0 ]
 *
 *  Reset law: x1^+ = rho1 * x1^-, x2^+ = rho2 * x2^-
 *  (potentially different reset ratios for each state)
 *
 *  Ref: Hazeleger, Heertjes, Nijmeijer (2011),
 *       "Second-order reset elements for stage control design"
 */
typedef struct {
    ResetSystem    *base;       /**< underlying reset system           */
    double          wn;         /**< natural frequency (rad/s)         */
    double          zeta;       /**< damping ratio                    */
    double          rho1;       /**< reset ratio for state 1           */
    double          rho2;       /**< reset ratio for state 2           */
    double          x0[2];      /**< initial state vector              */
    double          y;          /**< current output                    */
} SoreElement;

/* ================================================================
 * L1: Reset PID Controller (Zheng, Chait, Hollot, 2000)
 * ================================================================ */

/** ResetPID - PID controller with reset applied to the integrator.
 *
 *  The reset PID consists of:
 *    - Proportional term:  Kp * e(t)
 *    - Reset integrator:   KI * (Clegg or FORE integrator)
 *    - Derivative term:    Kd * (filtered derivative)
 *
 *  Transfer function form:
 *    C(s) = Kp + KI * R(s) + Kd * s / (tau_d * s + 1)
 *
 *  where R(s) is the describing function of the reset element
 *  (either Clegg integrator for full reset or FORE for partial reset).
 *
 *  The reset occurs on zero crossing of the error signal e(t),
 *  so the integrator is reset when the tracking error changes sign.
 *  This prevents integrator windup naturally without anti-windup
 *  schemes, while providing better phase characteristics.
 *
 *  Ref: [ZCH00]; Banos & Vidal (2005) on reset PI+CI compensator
 */
typedef struct {
    ResetSystem    *base_reset; /**< reset system for the I-part       */
    double          Kp;         /**< proportional gain                 */
    double          Ki;         /**< integral gain (multiplies reset output) */
    double          Kd;         /**< derivative gain                   */
    double          tau_d;      /**< derivative filter time const (N=1/tau_d) */
    double          i_state;    /**< integrator state (for output calc) */
    double          d_state;    /**< derivative filter state           */
    double          e_prev;     /**< previous error for derivative     */
    double          y;          /**< current PID output                */
    double          y_p;        /**< proportional part                 */
    double          y_i;        /**< integral (reset) part             */
    double          y_d;        /**< derivative part                   */
} ResetPID;

/* ================================================================
 * L1: Reset Lead-Lag Compensator
 * ================================================================ */

/** ResetLeadLag - Lead or lag compensator with reset capability.
 *
 *  Base transfer function (lead):
 *    G_lead(s) = K * (tau_z * s + 1) / (tau_p * s + 1),  tau_z > tau_p
 *
 *  Base transfer function (lag):
 *    G_lag(s)  = K * (tau_z * s + 1) / (tau_p * s + 1),  tau_z < tau_p
 *
 *  State-space realization (observable canonical form):
 *    A = -1/tau_p,  B = K*(tau_z/tau_p - 1)/tau_p,  C = 1,  D = K*tau_z/tau_p
 *
 *  Reset law: x^+ = rho * x^- (at zero crossing of input)
 *
 *  Ref: [BB12] Section 3.5
 */
typedef struct {
    ResetSystem    *base;       /**< underlying reset system           */
    double          K;          /**< gain                             */
    double          tau_z;      /**< zero time constant               */
    double          tau_p;      /**< pole time constant               */
    double          rho;        /**< reset ratio                      */
    double          x;          /**< compensator state                */
    double          y;          /**< current output                   */
    bool            is_lead;    /**< true if lead, false if lag       */
} ResetLeadLag;

/* ================================================================
 * L5: Element Lifecycle and Evaluation Functions
 * ================================================================ */

/* ---- Clegg Integrator ---- */

/** Create and initialize a Clegg integrator. */
CleggIntegrator* clegg_create(void);

/** Free a Clegg integrator. */
void clegg_free(CleggIntegrator *ci);

/** Evaluate the Clegg integrator for one step.
 *  dt: time step, e: current error, e_prev: previous error.
 *  Handles zero-crossing detection and reset internally.
 *  Returns the integrator output y. */
double clegg_step(CleggIntegrator *ci, double dt, double e, double e_prev);

/** Get the number of resets that have occurred. */
int clegg_num_resets(const CleggIntegrator *ci);

/** Compute the describing function magnitude at frequency w (rad/s).
 *  N_Clegg(w) = 1.62 / w (magnitude). */
double clegg_df_magnitude(double w);

/** Compute the describing function phase at frequency w.
 *  Returns phase in radians: -0.905 rad (-51.85 deg). */
double clegg_df_phase(void);

/* ---- FORE ---- */

/** Create a First Order Reset Element.
 *  K: DC gain, tau: time constant (>0), rho: reset ratio in [0,1).
 *  Returns NULL on invalid parameters. */
ForeElement* fore_create(double K, double tau, double rho);

/** Free a FORE element. */
void fore_free(ForeElement *fe);

/** Evaluate the FORE for one step.
 *  dt: time step, e: current error, e_prev: previous error.
 *  Handles zero-crossing detection and reset internally.
 *  Returns the filter output y. */
double fore_step(ForeElement *fe, double dt, double e, double e_prev);

/** Update the reset ratio. Changes the jump map accordingly.
 *  Returns false if rho is out of [0,1). */
bool fore_set_rho(ForeElement *fe, double rho);

/** Compute the describing function of FORE at frequency w.
 *  N_FORE(jw) = K / (jw*tau + 1) * (1 + jw*tau*rho) / (denominator...)
 *  Returns magnitude. */
double fore_df_magnitude(const ForeElement *fe, double w);

/** Compute the describing function phase of FORE at frequency w.
 *  Returns phase in radians. */
double fore_df_phase(const ForeElement *fe, double w);

/* ---- SORE ---- */

/** Create a Second Order Reset Element.
 *  wn: natural frequency (>0), zeta: damping ratio (>0),
 *  rho1, rho2: reset ratios for states 1 and 2. */
SoreElement* sore_create(double wn, double zeta, double rho1, double rho2);

/** Free a SORE element. */
void sore_free(SoreElement *se);

/** Evaluate the SORE for one step with ZC reset detection. */
double sore_step(SoreElement *se, double dt, double e, double e_prev);

/* ---- Reset PID ---- */

/** Create a reset PID controller.
 *  Kp, Ki, Kd: PID gains. tau_d: derivative filter constant.
 *  reset_rho: reset ratio for the integrator part (0 = Clegg). */
ResetPID* reset_pid_create(double Kp, double Ki, double Kd,
                            double tau_d, double reset_rho);

/** Free a reset PID controller. */
void reset_pid_free(ResetPID *pid);

/** Evaluate the reset PID for one step.
 *  dt: time step, e: current error, e_prev: previous error.
 *  Returns the total control output. */
double reset_pid_step(ResetPID *pid, double dt, double e, double e_prev);

/** Get the individual P, I (reset), and D components.
 *  Values placed in *yp, *yi, *yd respectively. */
void reset_pid_get_components(const ResetPID *pid, double *yp, double *yi, double *yd);

/** Reset the PID integrator manually (useful for initialization). */
void reset_pid_manual_reset(ResetPID *pid);

/* ---- Reset Lead-Lag ---- */

/** Create a reset lead-lag compensator.
 *  K: gain, tau_z: zero time constant, tau_p: pole time constant.
 *  rho: reset ratio. is_lead: true for lead, false for lag. */
ResetLeadLag* reset_leadlag_create(double K, double tau_z, double tau_p,
                                    double rho, bool is_lead);

/** Free a reset lead-lag compensator. */
void reset_leadlag_free(ResetLeadLag *rll);

/** Evaluate the reset lead-lag for one step. */
double reset_leadlag_step(ResetLeadLag *rll, double dt, double e, double e_prev);

/** Compute the frequency response (magnitude and phase) at frequency w.
 *  mag: output magnitude, phase: output phase in radians. */
void reset_leadlag_freqresp(const ResetLeadLag *rll, double w,
                             double *mag, double *phase);

#ifdef __cplusplus
}
#endif

#endif /* RESET_ELEMENT_H */