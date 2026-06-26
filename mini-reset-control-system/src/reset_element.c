/* reset_element.c - Reset Control Element Implementations
 *
 * Implements Clegg integrator, FORE, SORE, Reset PID, and
 * Reset Lead-Lag compensators as defined in reset_element.h.
 *
 * Each function implements an independent knowledge point:
 *   - clegg_step: L1/L5 Clegg integrator with ZC reset detection
 *   - clegg_df_magnitude: L2 describing function analysis
 *   - fore_step: L1/L5 FORE evaluation with partial reset
 *   - fore_df_magnitude: L2 FORE describing function
 *   - sore_step: L1 second-order reset element
 *   - reset_pid_step: L1/L5 reset PID controller
 *   - reset_leadlag_step: L1 reset lead-lag compensator
 *   - reset_leadlag_freqresp: L5 frequency response computation
 */

#include "reset_element.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ================================================================
 * L1: Clegg Integrator [Cle58]
 * ================================================================ */

CleggIntegrator* clegg_create(void)
{
    CleggIntegrator *ci = (CleggIntegrator*)malloc(sizeof(CleggIntegrator));
    if (!ci) return NULL;
    memset(ci, 0, sizeof(CleggIntegrator));

    ci->base = reset_clegg_create();
    if (!ci->base) {
        free(ci);
        return NULL;
    }
    ci->x0 = 0.0;
    ci->y  = 0.0;
    ci->n_resets = 0;
    return ci;
}

void clegg_free(CleggIntegrator *ci)
{
    if (!ci) return;
    reset_sys_free(ci->base);
    free(ci);
}

double clegg_step(CleggIntegrator *ci, double dt, double e, double e_prev)
{
    if (!ci || !ci->base) return 0.0;

    ResetSystem *rsys = ci->base;

    /* Check for zero-crossing and apply reset */
    ResetResult res = reset_check_and_apply(rsys, e, e_prev);
    if (res == RESET_OK) {
        ci->n_resets++;
    }

    /* Forward Euler integration: dx = e * dt */
    rsys->x_c[0] += e * dt;

    /* Update time and previous error */
    rsys->t += dt;
    rsys->e_prev = e;

    ci->y = rsys->x_c[0];
    return ci->y;
}

int clegg_num_resets(const CleggIntegrator *ci)
{
    if (!ci) return 0;
    return ci->n_resets;
}

/* Describing function of Clegg integrator:
 * N_Clegg(jw) = 1.62 / w * exp(-j * 51.85 deg)
 * Ref: [Cle58], [BB12] Eq. (3.4)
 */
double clegg_df_magnitude(double w)
{
    if (w <= 0.0) return INFINITY;
    return 1.62 / w;
}

double clegg_df_phase(void)
{
    /* -51.85 degrees = -0.905 radians */
    return -0.905;
}

/* ================================================================
 * L1: First Order Reset Element (FORE) [HR75]
 * ================================================================ */

ForeElement* fore_create(double K, double tau, double rho)
{
    if (tau <= 0.0 || rho < 0.0 || rho >= 1.0) return NULL;

    ForeElement *fe = (ForeElement*)malloc(sizeof(ForeElement));
    if (!fe) return NULL;
    memset(fe, 0, sizeof(ForeElement));

    fe->base = reset_fore_create(K, tau, rho);
    if (!fe->base) {
        free(fe);
        return NULL;
    }

    fe->K   = K;
    fe->tau = tau;
    fe->rho = rho;
    fe->x0  = 0.0;
    fe->y   = 0.0;

    return fe;
}

void fore_free(ForeElement *fe)
{
    if (!fe) return;
    reset_sys_free(fe->base);
    free(fe);
}

double fore_step(ForeElement *fe, double dt, double e, double e_prev)
{
    if (!fe || !fe->base) return 0.0;

    ResetSystem *rsys = fe->base;

    /* Check for zero-crossing and apply reset */
    (void)reset_check_and_apply(rsys, e, e_prev);
    /* Note: reset count tracked internally in rsys->interval */

    /* Forward Euler: dx = (-x/tau + K*e/tau) * dt */
    double x = rsys->x_c[0];
    double dx = (-x / fe->tau + fe->K * e / fe->tau);
    rsys->x_c[0] = x + dx * dt;

    rsys->t += dt;
    rsys->e_prev = e;

    fe->y = rsys->x_c[0];
    return fe->y;
}

bool fore_set_rho(ForeElement *fe, double rho)
{
    if (!fe || rho < 0.0 || rho >= 1.0) return false;
    fe->rho = rho;
    return (reset_set_ratio(fe->base, rho) == RESET_OK);
}

/* FORE Describing Function:
 * N_FORE(jw) = K/(jw*tau + 1) * correction_factor(rho, w)
 *
 * The correction factor accounts for the reset action and
 * depends on rho. For pure sinusoid input e = A*sin(wt),
 * the describing function provides the fundamental component.
 *
 * Ref: [BB12] Section 3.3.2, Eq. (3.25)-(3.28)
 */
double fore_df_magnitude(const ForeElement *fe, double w)
{
    if (!fe || w <= 0.0) return 0.0;

    double tau = fe->tau;
    double K   = fe->K;
    double rho = fe->rho;

    /* Base linear magnitude: |K / (jw*tau + 1)| */
    double denom = sqrt(1.0 + w * w * tau * tau);
    double mag_base = K / denom;

    /* Reset correction: approximate as (1+rho) factor
     * More precise formula involves the transcendental equation
     * from the harmonic balance analysis [BB12, Eq.3.26]
     */
    double correction = 2.0 / M_PI * (1.0 + rho);
    return mag_base * correction;
}

double fore_df_phase(const ForeElement *fe, double w)
{
    if (!fe || w <= 0.0) return 0.0;

    double tau = fe->tau;
    double rho = fe->rho;

    /* Base linear phase: -atan(w*tau) */
    double phase_base = -atan(w * tau);

    /* Reset phase lead: the reset provides extra phase advance
     * proportional to (1 - rho). For Clegg (rho=0): ~38.2 deg lead.
     * For FORE with rho>0: less phase lead. */
    double phase_lead = (1.0 - rho) * 0.667; /* ~38.2 deg max -> radians */

    return phase_base + phase_lead;
}

/* ================================================================
 * L1: Second Order Reset Element (SORE)
 *
 * Implementation of 2nd-order reset element with two independent
 * reset ratios for each state variable.
 *
 * Base TF: G(s) = wn^2 / (s^2 + 2*zeta*wn*s + wn^2)
 *
 * Ref: Hazeleger, Heertjes, Nijmeijer (2011)
 * ================================================================ */

SoreElement* sore_create(double wn, double zeta, double rho1, double rho2)
{
    if (wn <= 0.0 || zeta <= 0.0) return NULL;
    if (rho1 < 0.0 || rho1 >= 1.0) return NULL;
    if (rho2 < 0.0 || rho2 >= 1.0) return NULL;

    SoreElement *se = (SoreElement*)malloc(sizeof(SoreElement));
    if (!se) return NULL;
    memset(se, 0, sizeof(SoreElement));

    /* Create base reset system with 2 states, 1 input, 1 output */
    se->base = reset_sys_create(2);
    if (!se->base) { free(se); return NULL; }

    ResetLinearBase *flow = se->base->flow;
    flow->n = 2; flow->m = 1; flow->p = 1;

    /* Allocate flow matrices */
    flow->A = (double*)calloc(4, sizeof(double)); /* 2x2 */
    flow->B = (double*)calloc(2, sizeof(double)); /* 2x1 */
    flow->C = (double*)calloc(2, sizeof(double)); /* 1x2 */
    flow->D = (double*)calloc(1, sizeof(double)); /* 1x1 */

    if (!flow->A || !flow->B || !flow->C || !flow->D) {
        sore_free(se); return NULL;
    }

    /* Controllable canonical form:
     * A = [   0       1    ]
     *     [ -wn^2  -2*zeta*wn ]
     */
    flow->A[0 * 2 + 0] = 0.0;
    flow->A[0 * 2 + 1] = 1.0;
    flow->A[1 * 2 + 0] = -wn * wn;
    flow->A[1 * 2 + 1] = -2.0 * zeta * wn;

    /* B = [0; wn^2] */
    flow->B[0] = 0.0;
    flow->B[1] = wn * wn;

    /* C = [1, 0] (first state is output) */
    flow->C[0] = 1.0;
    flow->C[1] = 0.0;

    /* D = [0] */
    flow->D[0] = 0.0;

    /* Jump map: Ar = diag(rho1, rho2), Br = [0; 0] */
    se->base->jump->nc = 2;
    se->base->jump->Ar = (double*)calloc(4, sizeof(double));
    se->base->jump->Br = (double*)calloc(2, sizeof(double));

    if (!se->base->jump->Ar || !se->base->jump->Br) {
        sore_free(se); return NULL;
    }

    se->base->jump->Ar[0] = rho1;
    se->base->jump->Ar[3] = rho2;

    se->wn   = wn;
    se->zeta = zeta;
    se->rho1 = rho1;
    se->rho2 = rho2;
    se->x0[0] = 0.0;
    se->x0[1] = 0.0;
    se->y = 0.0;

    return se;
}

void sore_free(SoreElement *se)
{
    if (!se) return;
    reset_sys_free(se->base);
    free(se);
}

double sore_step(SoreElement *se, double dt, double e, double e_prev)
{
    if (!se || !se->base) return 0.0;

    ResetSystem *rsys = se->base;

    /* Zero-crossing detection and reset */
    reset_check_and_apply(rsys, e, e_prev);

    /* Euler integration for 2nd order system:
     * dx1 = x2 * dt
     * dx2 = (-wn^2*x1 - 2*zeta*wn*x2 + wn^2*e) * dt
     */
    double x1 = rsys->x_c[0];
    double x2 = rsys->x_c[1];
    double wn2 = se->wn * se->wn;

    rsys->x_c[0] = x1 + x2 * dt;
    rsys->x_c[1] = x2 + (-wn2 * x1 - 2.0 * se->zeta * se->wn * x2 + wn2 * e) * dt;

    rsys->t += dt;
    rsys->e_prev = e;

    se->y = rsys->x_c[0]; /* output = x1 */
    return se->y;
}

/* ================================================================
 * L1/L5: Reset PID Controller [ZCH00]
 *
 * PID with reset on the integrator (I) part only.
 * The P and D terms remain linear; only the I-term undergoes
 * reset on zero crossing of the error signal.
 *
 * This provides natural anti-windup and better phase
 * characteristics compared to linear PID.
 * ================================================================ */

ResetPID* reset_pid_create(double Kp, double Ki, double Kd,
                            double tau_d, double reset_rho)
{
    ResetPID *pid = (ResetPID*)malloc(sizeof(ResetPID));
    if (!pid) return NULL;
    memset(pid, 0, sizeof(ResetPID));

    /* Create FORE as the reset integrator:
     * For pure integrator behavior, use large tau and gain Ki.
     * The FORE K/(tau*s+1) approximates integrator for tau >> 1.
     * Alternatively for Clegg: use reset_clegg with Ki gain. */
    if (reset_rho < 1e-12) {
        /* Use Clegg integrator */
        pid->base_reset = reset_clegg_create();
    } else {
        /* Use FORE with large tau to approximate integrator */
        double tau_int = 1000.0;
        pid->base_reset = reset_fore_create(Ki, tau_int, reset_rho);
    }

    if (!pid->base_reset) { free(pid); return NULL; }

    pid->Kp     = Kp;
    pid->Ki     = Ki;
    pid->Kd     = Kd;
    pid->tau_d  = tau_d;
    pid->i_state = 0.0;
    pid->d_state = 0.0;
    pid->e_prev  = 0.0;
    pid->y       = 0.0;
    pid->y_p     = 0.0;
    pid->y_i     = 0.0;
    pid->y_d     = 0.0;

    return pid;
}

void reset_pid_free(ResetPID *pid)
{
    if (!pid) return;
    reset_sys_free(pid->base_reset);
    free(pid);
}

double reset_pid_step(ResetPID *pid, double dt, double e, double e_prev)
{
    if (!pid) return 0.0;

    /* --- Proportional term --- */
    pid->y_p = pid->Kp * e;

    /* --- Integral term (with reset) --- */
    if (pid->base_reset) {
        ResetSystem *rsys = pid->base_reset;

        /* For Clegg integrator: direct integration with reset
         * For FORE: use the fore dynamics */
        bool is_clegg = (fabs(rsys->ratio->rho) < 1e-12);
        if (is_clegg) {
            /* Check zero-crossing reset */
            (void)reset_check_and_apply(rsys, e, e_prev);
            /* Integrate: di_state = e * dt */
            rsys->x_c[0] += e * dt;
            rsys->t += dt;
            pid->i_state = rsys->x_c[0];
        } else {
            /* FORE integration */
            double x = rsys->x_c[0];
            double tau = 1000.0; /* should match creation */
            double dx = (-x / tau + pid->Ki * e / tau);
            (void)reset_check_and_apply(rsys, e, e_prev);
            rsys->x_c[0] = x + dx * dt;
            rsys->t += dt;
            pid->i_state = rsys->x_c[0];
        }
        rsys->e_prev = e;
    }
    pid->y_i = pid->i_state; /* Note: Ki already baked into the element */

    /* --- Derivative term (filtered) ---
     * D(s) = Kd * s / (tau_d * s + 1) applied to e
     * State-space: dx_d/dt = -(1/tau_d)*x_d + (Kd/tau_d)*e_dot
     * where e_dot ≈ (e - e_prev)/dt
     */
    if (pid->tau_d > 0.0 && dt > 0.0) {
        double e_dot = (e - pid->e_prev) / dt;
        double d_state = pid->d_state;
        double dd = (-d_state / pid->tau_d + pid->Kd * e_dot / pid->tau_d);
        pid->d_state = d_state + dd * dt;
        pid->y_d = pid->d_state;
    } else {
        pid->y_d = 0.0;
    }
    pid->e_prev = e;

    /* Total output */
    pid->y = pid->y_p + pid->y_i + pid->y_d;

    return pid->y;
}

void reset_pid_get_components(const ResetPID *pid, double *yp, double *yi, double *yd)
{
    if (!pid) return;
    if (yp) *yp = pid->y_p;
    if (yi) *yi = pid->y_i;
    if (yd) *yd = pid->y_d;
}

void reset_pid_manual_reset(ResetPID *pid)
{
    if (!pid || !pid->base_reset) return;
    ResetSystem *rsys = pid->base_reset;
    memset(rsys->x_c, 0, (size_t)rsys->nc * sizeof(double));
    pid->i_state = 0.0;
    pid->y_i = 0.0;
}

/* ================================================================
 * L1: Reset Lead-Lag Compensator
 *
 * Transforms a classical lead/lag network by adding reset capability
 * to the single state. This provides frequency-dependent phase
 * shaping with the added nonlinear phase advantage of reset.
 *
 * Base TF (observable canonical form):
 *   A = -1/tau_p
 *   B = K * (tau_z/tau_p - 1) / tau_p
 *   C = 1
 *   D = K * tau_z / tau_p
 * ================================================================ */

ResetLeadLag* reset_leadlag_create(double K, double tau_z, double tau_p,
                                    double rho, bool is_lead)
{
    if (tau_p <= 0.0 || tau_z <= 0.0) return NULL;
    if (rho < 0.0 || rho >= 1.0) return NULL;

    ResetLeadLag *rll = (ResetLeadLag*)malloc(sizeof(ResetLeadLag));
    if (!rll) return NULL;
    memset(rll, 0, sizeof(ResetLeadLag));

    rll->base = reset_sys_create(1);
    if (!rll->base) { free(rll); return NULL; }

    ResetLinearBase *flow = rll->base->flow;
    flow->n = 1; flow->m = 1; flow->p = 1;

    flow->A = (double*)calloc(1, sizeof(double));
    flow->B = (double*)calloc(1, sizeof(double));
    flow->C = (double*)calloc(1, sizeof(double));
    flow->D = (double*)calloc(1, sizeof(double));

    if (!flow->A || !flow->B || !flow->C || !flow->D) {
        reset_leadlag_free(rll); return NULL;
    }

    /* Observable canonical form:
     * A = -1/tau_p
     * B = K*(tau_z/tau_p - 1)/tau_p = K*(tau_z - tau_p)/(tau_p^2)
     * C = 1
     * D = K * tau_z / tau_p
     */
    flow->A[0] = -1.0 / tau_p;
    flow->B[0] = K * (tau_z / tau_p - 1.0) / tau_p;
    flow->C[0] = 1.0;
    flow->D[0] = K * tau_z / tau_p;

    /* Jump map: Ar = [rho], Br = [0] */
    rll->base->jump->nc = 1;
    rll->base->jump->Ar = (double*)calloc(1, sizeof(double));
    rll->base->jump->Br = (double*)calloc(1, sizeof(double));
    if (!rll->base->jump->Ar || !rll->base->jump->Br) {
        reset_leadlag_free(rll); return NULL;
    }
    rll->base->jump->Ar[0] = rho;

    rll->K       = K;
    rll->tau_z   = tau_z;
    rll->tau_p   = tau_p;
    rll->rho     = rho;
    rll->x       = 0.0;
    rll->y       = 0.0;
    rll->is_lead = is_lead;

    return rll;
}

void reset_leadlag_free(ResetLeadLag *rll)
{
    if (!rll) return;
    reset_sys_free(rll->base);
    free(rll);
}

double reset_leadlag_step(ResetLeadLag *rll, double dt, double e, double e_prev)
{
    if (!rll || !rll->base) return 0.0;

    ResetSystem *rsys = rll->base;

    /* Zero-crossing detection and reset */
    reset_check_and_apply(rsys, e, e_prev);

    /* Forward Euler: dx = (A*x + B*e) * dt */
    double x = rsys->x_c[0];
    double dx = rsys->flow->A[0] * x + rsys->flow->B[0] * e;
    rsys->x_c[0] = x + dx * dt;

    /* Output: y = C*x + D*e */
    rll->x = rsys->x_c[0];
    rll->y = rsys->flow->C[0] * rll->x + rsys->flow->D[0] * e;

    rsys->t += dt;
    rsys->e_prev = e;

    return rll->y;
}

void reset_leadlag_freqresp(const ResetLeadLag *rll, double w,
                             double *mag, double *phase)
{
    if (!rll || !mag || !phase) return;

    double K     = rll->K;
    double tau_z = rll->tau_z;
    double tau_p = rll->tau_p;

    /* G(jw) = K * (jw*tau_z + 1) / (jw*tau_p + 1) */
    double num_real = 1.0;
    double num_imag = w * tau_z;
    double den_real = 1.0;
    double den_imag = w * tau_p;

    double num_mag = sqrt(num_real * num_real + num_imag * num_imag);
    double den_mag = sqrt(den_real * den_real + den_imag * den_imag);

    *mag   = K * num_mag / den_mag;
    *phase = atan2(num_imag, num_real) - atan2(den_imag, den_real);
}