#ifndef DTA_STATE_DEPENDENT_H
#define DTA_STATE_DEPENDENT_H
#include "dta_core.h"
#include "dta_switch_signal.h"

/* ==============================================================
 * dta_state_dependent.h - State-Dependent Switching
 *
 * State-dependent switching defines σ(x(t)) based on the current
 * state, rather than σ(t) as a predetermined function of time.
 *
 * Common approach: partition state space R^n into regions Ω_i,
 * and use mode i when x ∈ Ω_i.
 *
 * Hysteresis switching (Morse 1996):
 *   Add hysteresis bands to prevent chattering (infinite fast
 *   switching) at region boundaries. A dwell-time naturally
 *   emerges from the hysteresis width and state velocity.
 *
 * Key theorem: If there exists a family of Lyapunov functions
 *   {V_i} such that V̇_i(x) < 0 on Ω_i and V_j(x) ≤ V_i(x) at
 *   switching boundaries, then the system is GAS.
 *
 * References:
 *   Morse (1996) IEEE TAC 41(10):1413-1431
 *   Liberzon (2003) Ch. 3.4 — State-dependent switching
 *   Pettersson (2003) "Synthesis of switched linear systems",
 *     IEEE CDC
 * ============================================================== */

/* --- State space partitioning --- */
typedef enum {
    DTA_PARTITION_QUADRANT = 0,
    DTA_PARTITION_HALFPLANE = 1,
    DTA_PARTITION_POLYTOPE = 2,
    DTA_PARTITION_VORONOI = 3,
    DTA_PARTITION_ELLIPSOID = 4
} DTA_PartitionType;

/* --- A region in state space --- */
typedef struct {
    int region_id;
    DTA_PartitionType type;
    double* normal;            /* For halfplane: {x : n^T x ≤ b} */
    double bias;
    double* A_poly;            /* For polytope: A x ≤ b, m×n */
    double* b_poly;
    int m_constraints;         /* Number of halfplane constraints */
    double* center;            /* Center (for Voronoi/ellipsoid) */
    double* Q_ellip;           /* Shape matrix for ellipsoid: (x-c)^T Q (x-c) ≤ 1 */
    int n;                     /* State dimension */
} DTA_StateRegion;

/* --- State-dependent switching law --- */
typedef struct {
    DTA_StateRegion* regions;
    int n_regions;
    int* region_to_mode;       /* Map region → mode index */
    int n;                     /* State dimension */
    double hysteresis_band;    /* Width of hysteresis zone */
    double min_dwell;          /* Inherent minimum dwell time */
    bool use_hysteresis;
} DTA_StateSwitchLaw;

/* --- Hysteresis switching monitor --- */
typedef struct {
    int current_mode;
    int previous_mode;
    double last_switch_time;
    double hysteresis_time;    /* Time spent in hysteresis zone */
    bool in_hysteresis;
    bool chattering_detected;
} DTA_HysteresisMonitor;

/* --- API --- */

/** Create a state-dependent switching law */
DTA_StateSwitchLaw* dta_sdsl_create(int n_regions, int n,
                                     double hysteresis_band);
void dta_sdsl_free(DTA_StateSwitchLaw* law);

/** Define a halfplane region: {x : n^T x ≤ b} */
int dta_sdsl_set_halfplane_region(DTA_StateSwitchLaw* law,
    int region_id, const double* normal, double bias, int mode);

/** Define a polytope region: {x : A x ≤ b} */
int dta_sdsl_set_polytope_region(DTA_StateSwitchLaw* law,
    int region_id, const double* A_poly, const double* b_poly,
    int m_constraints, int mode);

/** Define a quadrant region based on sign of coordinates */
int dta_sdsl_set_quadrant_region(DTA_StateSwitchLaw* law,
    int region_id, const int* signs, int mode);

/** Determine the active mode for state x based on the partition.
 *  With hysteresis, returns previous mode if x is in hysteresis zone. */
int dta_sdsl_active_mode(const DTA_StateSwitchLaw* law,
                          const double* x, double t, int prev_mode);

/** Check if state x is in region r */
bool dta_sdsl_in_region(const DTA_StateRegion* region, const double* x);

/** Compute distance from state x to region boundary */
double dta_sdsl_boundary_distance(const DTA_StateRegion* region,
                                   const double* x);

/** Simulate state-dependent switching from initial x0.
 *  Applies the switching law at each step. */
DTA_StateTrajectory* dta_sdsl_simulate(const DTA_SwitchedSystem* sys,
    const DTA_StateSwitchLaw* law, const double* x0,
    double t_end, double dt);

/** Detect chattering: whether switches occur too frequently.
 *  Returns true if the minimum inter-switch time < min_dwell. */
bool dta_sdsl_detect_chattering(const DTA_SwitchingSignal* sig,
                                 double min_dwell);

/** Design hysteresis to prevent chattering for given system.
 *  Computes minimum hysteresis band width to guarantee
 *  minimum dwell time τ_d. */
double dta_sdsl_design_hysteresis(const DTA_SwitchedSystem* sys,
                                   const DTA_StateSwitchLaw* law,
                                   double tau_d);

/** Create hysteresis monitor */
DTA_HysteresisMonitor dta_sdsl_monitor_create(void);

/** Update hysteresis monitor with new state information */
void dta_sdsl_monitor_update(DTA_HysteresisMonitor* mon,
    const double* x, double t, int current_mode);

/** Voronoi partition: assign mode based on nearest center */
int dta_sdsl_voronoi_mode(const double* centers, int n_centers,
                           int n, const double* x);

/** Ellipsoidal guard condition: switch when (x-c)^T Q (x-c) = 1 */
bool dta_sdsl_guard_ellipsoid(const double* x, const double* center,
                               const double* Q, int n, double tol);

/** Scale-independent hysteresis: dwell time scales with state norm */
double dta_sdsl_scale_independent_hysteresis(const DTA_SwitchedSystem* sys,
    int from_mode, int to_mode, double state_norm);

#endif /* DTA_STATE_DEPENDENT_H */
