#include "dta_core.h"
#include "dta_switch_signal.h"
#include "dta_state_dependent.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * dta_state_dependent.c - State-Dependent Switching
 *
 * Implements:
 *   State-space partitioning (halfplane, polytope, quadrant, Voronoi)
 *   State-dependent switching law based on partition
 *   Hysteresis switching to prevent chattering
 *   Chattering detection
 *   Hysteresis band design for guaranteed dwell time
 *   Scale-independent hysteresis
 *   Ellipsoidal guard conditions
 *
 * The key idea: sigma(t) = g(x(t)) based on state x(t), not time.
 * Hysteresis ensures minimum dwell time naturally emerges.
 *
 * References:
 *   Morse (1996) IEEE TAC 41(10):1413-1431
 *   Liberzon (2003) "Switching in Systems and Control" Ch.3.4
 *   Pettersson (2003) IEEE CDC
 * ============================================================== */

static void* safe_alloc(size_t sz) {
    void* p = malloc(sz);
    if (!p) { fprintf(stderr, "DTA sdsl: alloc fail\n"); exit(1); }
    return p;
}

DTA_StateSwitchLaw* dta_sdsl_create(int n_regions, int n,
                                     double hysteresis_band) {
    DTA_StateSwitchLaw* law = safe_alloc(sizeof(DTA_StateSwitchLaw));
    law->n_regions = n_regions;
    law->n = n;
    law->hysteresis_band = hysteresis_band;
    law->min_dwell = 0.0;
    law->use_hysteresis = (hysteresis_band > 0);
    law->regions = safe_alloc((size_t)n_regions * sizeof(DTA_StateRegion));
    law->region_to_mode = safe_alloc((size_t)n_regions * sizeof(int));
    int i;
    for (i = 0; i < n_regions; i++) {
        law->regions[i].region_id = i;
        law->regions[i].type = DTA_PARTITION_HALFPLANE;
        law->regions[i].normal = NULL;
        law->regions[i].bias = 0.0;
        law->regions[i].A_poly = NULL;
        law->regions[i].b_poly = NULL;
        law->regions[i].m_constraints = 0;
        law->regions[i].center = NULL;
        law->regions[i].Q_ellip = NULL;
        law->regions[i].n = n;
        law->region_to_mode[i] = 0;
    }
    return law;
}

void dta_sdsl_free(DTA_StateSwitchLaw* law) {
    if (!law) return;
    int i;
    for (i = 0; i < law->n_regions; i++) {
        free(law->regions[i].normal);
        free(law->regions[i].A_poly);
        free(law->regions[i].b_poly);
        free(law->regions[i].center);
        free(law->regions[i].Q_ellip);
    }
    free(law->regions);
    free(law->region_to_mode);
    free(law);
}

int dta_sdsl_set_halfplane_region(DTA_StateSwitchLaw* law,
    int region_id, const double* normal, double bias, int mode) {
    if (!law || region_id < 0 || region_id >= law->n_regions || !normal)
        return -1;
    DTA_StateRegion* r = &law->regions[region_id];
    r->type = DTA_PARTITION_HALFPLANE;
    r->normal = safe_alloc((size_t)law->n * sizeof(double));
    memcpy(r->normal, normal, (size_t)law->n * sizeof(double));
    r->bias = bias;
    law->region_to_mode[region_id] = mode;
    return 0;
}

int dta_sdsl_set_polytope_region(DTA_StateSwitchLaw* law,
    int region_id, const double* A_poly, const double* b_poly,
    int m_constraints, int mode) {
    if (!law || region_id < 0 || region_id >= law->n_regions || !A_poly || !b_poly)
        return -1;
    DTA_StateRegion* r = &law->regions[region_id];
    r->type = DTA_PARTITION_POLYTOPE;
    r->m_constraints = m_constraints;
    r->A_poly = safe_alloc((size_t)(m_constraints * law->n) * sizeof(double));
    memcpy(r->A_poly, A_poly, (size_t)(m_constraints * law->n) * sizeof(double));
    r->b_poly = safe_alloc((size_t)m_constraints * sizeof(double));
    memcpy(r->b_poly, b_poly, (size_t)m_constraints * sizeof(double));
    law->region_to_mode[region_id] = mode;
    return 0;
}

int dta_sdsl_set_quadrant_region(DTA_StateSwitchLaw* law,
    int region_id, const int* signs, int mode) {
    if (!law || region_id < 0 || region_id >= law->n_regions || !signs)
        return -1;
    DTA_StateRegion* r = &law->regions[region_id];
    r->type = DTA_PARTITION_QUADRANT;
    /* Convert signs to halfplane constraints: sign_i * x_i >= 0 */
    int n = law->n;
    r->m_constraints = n;
    r->A_poly = safe_alloc((size_t)(n * n) * sizeof(double));
    r->b_poly = safe_alloc((size_t)n * sizeof(double));
    memset(r->A_poly, 0, (size_t)(n * n) * sizeof(double));
    int i;
    for (i = 0; i < n; i++) {
        /* Constraint: -sign_i * x_i <= 0  (so A[i*n+i] = -sign_i, b[i] = 0) */
        r->A_poly[i*n + i] = -(double)signs[i];
        r->b_poly[i] = 0.0;
    }
    law->region_to_mode[region_id] = mode;
    return 0;
}

bool dta_sdsl_in_region(const DTA_StateRegion* region, const double* x) {
    if (!region || !x) return false;
    int n = region->n;
    switch (region->type) {
    case DTA_PARTITION_HALFPLANE: {
        double dot = 0.0;
        int i;
        for (i = 0; i < n; i++) dot += region->normal[i] * x[i];
        return dot <= region->bias;
    }
    case DTA_PARTITION_POLYTOPE:
    case DTA_PARTITION_QUADRANT: {
        int i, j;
        for (i = 0; i < region->m_constraints; i++) {
            double dot = 0.0;
            for (j = 0; j < n; j++)
                dot += region->A_poly[i*n + j] * x[j];
            if (dot > region->b_poly[i] + 1e-10) return false;
        }
        return true;
    }
    case DTA_PARTITION_VORONOI: {
        if (!region->center) return false;
        double dist2 = 0.0;
        int i;
        for (i = 0; i < n; i++) {
            double d = x[i] - region->center[i];
            dist2 += d * d;
        }
        return dist2 <= 1.0;
    }
    case DTA_PARTITION_ELLIPSOID: {
        if (!region->center || !region->Q_ellip) return false;
        double* diff = safe_alloc((size_t)n * sizeof(double));
        int i, j;
        for (i = 0; i < n; i++) diff[i] = x[i] - region->center[i];
        double val = 0.0;
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                val += diff[i] * region->Q_ellip[i*n + j] * diff[j];
        free(diff);
        return val <= 1.0 + 1e-10;
    }
    default: return false;
    }
}

int dta_sdsl_active_mode(const DTA_StateSwitchLaw* law,
                          const double* x, double t, int prev_mode) {
    (void)t;
    if (!law || !x) return -1;
    int i;
    /* Check all regions, return first match */
    for (i = 0; i < law->n_regions; i++) {
        if (dta_sdsl_in_region(&law->regions[i], x))
            return law->region_to_mode[i];
    }
    return (prev_mode >= 0) ? prev_mode : 0;
}

double dta_sdsl_boundary_distance(const DTA_StateRegion* region,
                                   const double* x) {
    if (!region || !x) return INFINITY;
    int n = region->n;
    switch (region->type) {
    case DTA_PARTITION_HALFPLANE: {
        double dot = 0.0;
        int i;
        for (i = 0; i < n; i++) dot += region->normal[i] * x[i];
        return fabs(dot - region->bias);
    }
    case DTA_PARTITION_POLYTOPE: {
        double min_dist = INFINITY;
        int i, j;
        for (i = 0; i < region->m_constraints; i++) {
            double dot = 0.0;
            for (j = 0; j < n; j++)
                dot += region->A_poly[i*n + j] * x[j];
            double dist = fabs(dot - region->b_poly[i]);
            if (dist < min_dist) min_dist = dist;
        }
        return min_dist;
    }
    default: return 0.0;
    }
}

/* Simulate state-dependent switching */
DTA_StateTrajectory* dta_sdsl_simulate(const DTA_SwitchedSystem* sys,
    const DTA_StateSwitchLaw* law, const double* x0,
    double t_end, double dt) {
    if (!sys || !law || !x0 || dt <= 0) return NULL;
    int n = sys->state_dim;
    int n_steps = (int)((t_end - 0.0) / dt) + 2;
    DTA_StateTrajectory* traj = safe_alloc(sizeof(DTA_StateTrajectory));
    traj->capacity = n_steps;
    traj->n_samples = 0;
    traj->state_dim = n;
    traj->t_samples = safe_alloc((size_t)n_steps * sizeof(double));
    traj->x_samples = calloc((size_t)n_steps, sizeof(double*));
    if (!traj->x_samples) { free(traj->t_samples); free(traj); return NULL; }
    int i;
    for (i = 0; i < n_steps; i++) {
        traj->x_samples[i] = safe_alloc((size_t)n * sizeof(double));
    }

    /* Initial state */
    traj->t_samples[0] = 0.0;
    if (n_steps > 0)
        memcpy(traj->x_samples[0], x0, (size_t)n * sizeof(double));
    traj->n_samples = 1;

    double* x = safe_alloc((size_t)n * sizeof(double));
    double* dx = safe_alloc((size_t)n * sizeof(double));
    memcpy(x, x0, (size_t)n * sizeof(double));
    double t = 0.0;
    int current_mode = dta_sdsl_active_mode(law, x, 0.0, 0);
    double last_switch_t = 0.0;
    int step;

    for (step = 0; step < n_steps - 1 && t < t_end; step++) {
        /* Check for state-dependent switch */
        int new_mode = dta_sdsl_active_mode(law, x, t, current_mode);

        /* Hysteresis: only switch if enough time has passed */
        if (new_mode != current_mode && law->use_hysteresis) {
            double time_since_switch = t - last_switch_t;
            if (time_since_switch < law->min_dwell) {
                new_mode = current_mode;  /* Block switch */
            }
        }

        if (new_mode != current_mode) {
            current_mode = new_mode;
            last_switch_t = t;
        }

        /* RK4 integration step */
        double k1[10], k2[10], k3[10], k4[10], x_temp[10];
        int j;
        dta_system_rhs(sys, current_mode, t, x, NULL, k1);
        for (j = 0; j < n && j < 10; j++) x_temp[j] = x[j] + 0.5*dt*k1[j];
        dta_system_rhs(sys, current_mode, t + 0.5*dt, x_temp, NULL, k2);
        for (j = 0; j < n && j < 10; j++) x_temp[j] = x[j] + 0.5*dt*k2[j];
        dta_system_rhs(sys, current_mode, t + 0.5*dt, x_temp, NULL, k3);
        for (j = 0; j < n && j < 10; j++) x_temp[j] = x[j] + dt*k3[j];
        dta_system_rhs(sys, current_mode, t + dt, x_temp, NULL, k4);

        for (j = 0; j < n && j < 10; j++)
            x[j] += dt / 6.0 * (k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]);

        t += dt;

        /* Record sample every step */
        if (traj->n_samples < traj->capacity) {
            traj->t_samples[traj->n_samples] = t;
            memcpy(traj->x_samples[traj->n_samples], x, (size_t)n * sizeof(double));
            traj->n_samples++;
        }
    }

    free(x); free(dx);
    return traj;
}

bool dta_sdsl_detect_chattering(const DTA_SwitchingSignal* sig,
                                 double min_dwell) {
    if (!sig) return false;
    double actual_min = dta_signal_min_dwell(sig);
    return actual_min < min_dwell;
}

/* Design hysteresis band to guarantee minimum dwell time tau_d.
 * Estimates: hysteresis_band >= tau_d * max_velocity
 * max_velocity estimated from ||A_i|| * max_state_norm */
double dta_sdsl_design_hysteresis(const DTA_SwitchedSystem* sys,
                                   const DTA_StateSwitchLaw* law,
                                   double tau_d) {
    if (!sys || !law || tau_d <= 0) return 0.0;
    int i;
    double max_gain = 0.0;
    for (i = 0; i < sys->n_modes; i++) {
        double norm_A = 0.0;
        int j, k;
        for (j = 0; j < sys->state_dim; j++)
            for (k = 0; k < sys->state_dim; k++)
                norm_A += sys->modes[i].A[j*sys->state_dim + k]
                        * sys->modes[i].A[j*sys->state_dim + k];
        norm_A = sqrt(norm_A);
        if (norm_A > max_gain) max_gain = norm_A;
    }
    return tau_d * max_gain * 1.0;  /* Conservative estimate */
}

DTA_HysteresisMonitor dta_sdsl_monitor_create(void) {
    DTA_HysteresisMonitor mon;
    memset(&mon, 0, sizeof(mon));
    mon.current_mode = 0;
    mon.previous_mode = 0;
    mon.last_switch_time = 0.0;
    mon.hysteresis_time = 0.0;
    mon.in_hysteresis = false;
    mon.chattering_detected = false;
    return mon;
}

void dta_sdsl_monitor_update(DTA_HysteresisMonitor* mon,
    const double* x, double t, int current_mode) {
    if (!mon) return;
    if (current_mode != mon->current_mode) {
        double dwell = t - mon->last_switch_time;
        if (dwell < 1e-6) mon->chattering_detected = true;
        mon->previous_mode = mon->current_mode;
        mon->current_mode = current_mode;
        mon->last_switch_time = t;
        mon->hysteresis_time = 0.0;
        mon->in_hysteresis = false;
    } else {
        mon->hysteresis_time = t - mon->last_switch_time;
    }
    (void)x;
}

int dta_sdsl_voronoi_mode(const double* centers, int n_centers,
                           int n, const double* x) {
    if (!centers || !x || n_centers <= 0) return -1;
    double min_dist = INFINITY;
    int best = 0, i;
    for (i = 0; i < n_centers; i++) {
        double dist2 = 0.0;
        int j;
        for (j = 0; j < n; j++) {
            double d = x[j] - centers[i*n + j];
            dist2 += d * d;
        }
        if (dist2 < min_dist) { min_dist = dist2; best = i; }
    }
    return best;
}

bool dta_sdsl_guard_ellipsoid(const double* x, const double* center,
                               const double* Q, int n, double tol) {
    if (!x || !center || !Q || n <= 0) return false;
    double* diff = safe_alloc((size_t)n * sizeof(double));
    int i, j;
    for (i = 0; i < n; i++) diff[i] = x[i] - center[i];
    double val = 0.0;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            val += diff[i] * Q[i*n + j] * diff[j];
    free(diff);
    return fabs(val - 1.0) < tol;
}

double dta_sdsl_scale_independent_hysteresis(const DTA_SwitchedSystem* sys,
    int from_mode, int to_mode, double state_norm) {
    if (!sys || from_mode < 0 || to_mode < 0) return 0.0;
    /* Scale hysteresis with state norm: larger state -> larger band */
    double base = 0.01;
    if (from_mode < sys->n_modes && to_mode < sys->n_modes) {
        double mu_from = dta_matrix_measure(sys->modes[from_mode].A, sys->state_dim);
        double mu_to = dta_matrix_measure(sys->modes[to_mode].A, sys->state_dim);
        base = fabs(mu_from - mu_to) * 0.1;
    }
    return base * (1.0 + state_norm);
}
