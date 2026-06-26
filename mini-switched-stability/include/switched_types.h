#ifndef SWITCHED_TYPES_H
#define SWITCHED_TYPES_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============================================================================
 * Switched System Core Types (L1-L3)
 *
 * Based on foundational works in switched and hybrid systems:
 *   Liberzon, D. (2003). Switching in Systems and Control. Birkhauser.
 *   Branicky, M.S. (1998). Multiple Lyapunov functions and other analysis
 *     tools for switched and hybrid systems. IEEE TAC, 43(4), 475-482.
 *   Hespanha, J.P. & Morse, A.S. (1999). Stability of switched systems
 *     with average dwell-time. IEEE CDC.
 *   Morse, A.S. (1996). Supervisory control of families of linear
 *     set-point controllers. IEEE TAC, 41(10), 1413-1431.
 *   Sun, Z. & Ge, S.S. (2011). Stability Theory of Switched Dynamical
 *     Systems. Springer.
 * ============================================================================ */

/* --------------------------------------------------------------------------
 * L1: Core Definitions - Switched System Components
 * -------------------------------------------------------------------------- */

/** Subsystem mode index. Finite set of modes P = {0, 1, ..., p-1}. */
typedef int SubsystemMode;

/**
 * Switching signal type enumeration.
 * Describes the nature of the switching law sigma(t).
 */
typedef enum {
    SSIG_ARBITRARY    = 0,  /* sigma(t) can switch at any time, no constraints   */
    SSIG_PERIODIC     = 1,  /* sigma(t) follows a fixed periodic schedule         */
    SSIG_STATE_DEP    = 2,  /* sigma(t) = phi(x(t)), depending on state           */
    SSIG_TIME_DEP     = 3,  /* sigma(t) = psi(t), depending on time only          */
    SSIG_EVENT_TRIG   = 4,  /* sigma(t) switches on discrete events               */
    SSIG_MARKOVIAN    = 5,  /* sigma(t) is a finite-state Markov process          */
    SSIG_DWELL_TIME   = 6,  /* sigma(t) constrained by minimum dwell time tau_d   */
    SSIG_AVG_DWELL    = 7   /* sigma(t) constrained by average dwell time tau_a   */
} SwitchingSignalType;

/**
 * Stability mode enumeration for switched systems.
 */
typedef enum {
    SSTAB_UNSTABLE              = 0,  /* Not stable under given switching         */
    SSTAB_LYAPUNOV_STABLE       = 1,  /* Lyapunov stable (eps-delta definition)    */
    SSTAB_ASYMPTOTIC_STABLE     = 2,  /* Asymptotically stable                     */
    SSTAB_EXPONENTIAL_STABLE    = 3,  /* ||x(t)|| <= c e^{-lambda t} ||x(0)||      */
    SSTAB_GUES                  = 4,  /* Globally Uniformly Exponentially Stable   */
    SSTAB_PRACTICAL_STABLE      = 5,  /* Stable within a bounded region            */
    SSTAB_INPUT_TO_STATE        = 6   /* ISS: stability with respect to inputs     */
} SwitchedStabilityType;

/**
 * Lyapunov function type for a switched system.
 */
typedef enum {
    CLF_COMMON           = 0,  /* Common Lyapunov Function (CLF): V(x) = x^T P x   */
    CLF_MULTIPLE         = 1,  /* Multiple Lyapunov Functions: V_i(x) = x^T P_i x  */
    CLF_PIECEWISE        = 2,  /* Piecewise Lyapunov Function                       */
    CLF_DWELL_TIME       = 3,  /* Dwell-time Lyapunov: V_i decreases by factor mu   */
    CLF_AVERAGE_DWELL    = 4,  /* Average dwell-time Lyapunov: averaged decay       */
    CLF_PATH_COMPLETE    = 5   /* Path-complete Lyapunov functions                  */
} LyapunovFunctionType;

/* --------------------------------------------------------------------------
 * L2: Core Concepts - Vector and Matrix Types
 * -------------------------------------------------------------------------- */

/** Real vector (dense, contiguous storage). State x in R^n. */
typedef struct {
    double *data;
    int     n;
} SwitchedVector;

/** Real matrix (dense, row-major storage). A_i in R^{n x n}. */
typedef struct {
    double *data;
    int     rows;
    int     cols;
} SwitchedMatrix;

/**
 * A single linear subsystem of the switched system.
 * Dynamics: dx/dt = A_k x + B_k u.
 */
typedef struct {
    int             mode_id;
    SwitchedMatrix  A;
    SwitchedMatrix  B;
    double          spectral_radius;
    bool            is_hurwitz;
    bool            is_schur;
    char            description[64];
} SwitchedSubsystem;

/**
 * Switching signal sigma(t).
 * Describes which subsystem is active at time t.
 */
typedef struct {
    SwitchingSignalType type;
    int                *mode_sequence;
    double             *switch_times;
    int                 n_switches;
    int                 capacity;
    int                 current_mode;
    double              total_time;
} SwitchingSignal;

/* --------------------------------------------------------------------------
 * L3: Mathematical Structures - Switched System Model
 * -------------------------------------------------------------------------- */

/**
 * Complete switched system model.
 *
 * Dynamics: dx/dt(t) = f_{sigma(t)}(x(t), t)
 * For linear subsystems: dx/dt(t) = A_{sigma(t)} x(t)
 *
 * The switching signal sigma: [0, inf) -> P selects which subsystem
 * dynamics govern the evolution at each instant.
 */
typedef struct {
    char                name[128];
    int                 state_dim;
    int                 n_modes;
    SwitchedSubsystem **modes;
    SwitchingSignal    *signal;
    SwitchedVector      state;
    double              t;

    /* Stability analysis data */
    SwitchedStabilityType stability;
    double              *eigenvalues;
    int                  eigenvalue_dim;

    /* Lyapunov function data */
    LyapunovFunctionType  lf_type;
    SwitchedMatrix      **lyap_matrices;
    double               *lyap_values;
    bool                 *mode_active;

    /* Dwell time tracking */
    double   avg_dwell_time;
    double   min_dwell_time;
    int     *dwell_samples;
    int      dwell_count;
} SwitchedSystem;

/* --------------------------------------------------------------------------
 * L2: Lyapunov Function Structures
 * -------------------------------------------------------------------------- */

/** Common Lyapunov Function (CLF). V(x) = x^T P x, P = P^T > 0. */
typedef struct {
    SwitchedMatrix P;
    double         min_eig;
    double         max_eig;
    bool           is_valid;
    double         decay_rate;
} CommonLyapunovFunction;

/** Multiple Lyapunov Functions (MLF). Each mode i has V_i(x) = x^T P_i x. */
typedef struct {
    SwitchedMatrix  *P;
    double          *min_eig;
    double          *max_eig;
    bool            *is_valid;
    double          *decay_rates;
    int              n_modes;
    double           mu;
    bool             mlf_condition;
} MultipleLyapunovFunctions;

/**
 * Dwell-time based Lyapunov analysis.
 * Under minimum dwell time tau_d > 0, stability is guaranteed
 * if each subsystem is exponentially stable and switches are
 * sufficiently slow (separated by at least tau_d).
 */
typedef struct {
    double          tau_d;
    double          tau_a;
    int             N0;
    bool            slow_enough;
    double          stability_margin;
    double          mu_lyap;
    double          required_tau_a;
} DwellTimeAnalysis;

/* --------------------------------------------------------------------------
 * L3: Lie-Algebraic Stability Condition Types
 * -------------------------------------------------------------------------- */

/**
 * Lie-algebraic condition result.
 * If the Lie algebra generated by {A_0,...,A_{p-1}} is solvable,
 * then a CLF exists and GUES is guaranteed under arbitrary switching.
 *
 * Reference: Liberzon, Hespanha & Morse (1999). Stability of switched
 *   systems: a Lie-algebraic condition. Systems & Control Letters.
 */
typedef struct {
    bool     is_solvable;
    bool     is_nilpotent;
    bool     pair_commute;
    bool     simultaneously_triang;
    double   *lie_brackets;
    int      n_modes;
    int      dim;
} LieAlgebraCondition;

/* --------------------------------------------------------------------------
 * L6: Canonical Problem Types - Application Domains
 * -------------------------------------------------------------------------- */

/**
 * DC-DC converter switching topology (power electronics).
 * Two modes: switch ON (charging inductor) / switch OFF (discharging).
 * States: iL = inductor current, vC = capacitor voltage.
 *
 * ON mode:  diL/dt = (Vin - vC)/L,     dvC/dt = (iL - vC/R)/C
 * OFF mode: diL/dt = -vC/L,             dvC/dt = (iL - vC/R)/C
 */
typedef struct {
    double Vin;
    double Vout;
    double L;
    double C;
    double R_load;
    double duty_cycle;
    double freq;
    bool   mode_on;
    double iL;
    double vC;
} DCDCConverter;

/**
 * Thermostat control system (bang-bang control with hysteresis).
 * State: temperature and rate of change.
 * Mode switches: OFF, HEATING, COOLING based on deadband.
 */
typedef struct {
    double temp;
    double setpoint;
    double deadband;
    double heating_rate;
    double cooling_rate;
    double ambient_loss;
    int    current_mode;
    double t_since_switch;
} ThermostatSystem;

/**
 * Automated highway vehicle spacing control.
 * Mode switches as lead vehicle changes speed.
 */
typedef struct {
    double ego_speed;
    double lead_speed;
    double gap;
    double safe_gap;
    int    speed_mode;
    double control_input;
} VehicleSpacingControl;

/**
 * Networked control system with packet dropouts.
 * Dropouts cause mode switches between controller structures.
 */
typedef struct {
    double      *state;
    int          state_dim;
    double      *ctrl_gain;
    bool         packet_received;
    double       dropout_rate;
    int          consecutive_losses;
    double       max_allowable_loss;
    bool         is_stable;
} NetworkedControlDropout;

/* --------------------------------------------------------------------------
 * L3: Numerical solver types
 * -------------------------------------------------------------------------- */

/** ODE solver configuration for switched systems. */
typedef struct {
    double  dt;
    double  t_end;
    int     max_steps;
    double  event_tol;
    bool    detect_switches;
} SolverConfig;

/* --------------------------------------------------------------------------
 * L5: Algorithm structure types
 * -------------------------------------------------------------------------- */

/**
 * Linear Matrix Inequality (LMI) feasibility problem.
 * Form: F(x) = F_0 + sum_i x_i F_i > 0
 * Used for Lyapunov matrix computation.
 */
typedef struct {
    SwitchedMatrix  *F;
    int              n_terms;
    int              dim;
    double          *solution;
    bool             found;
    int              iters;
} LMIProblem;

/** Eigenvalue result for stability analysis. */
typedef struct {
    double  real;
    double  imag;
    double  magnitude;
    bool    stable;
} EigenvalueResult;

#define QR_MAX_ITERS 1000
#define QR_TOL 1e-12

/** QR algorithm workspace for eigenvalue computation. */
typedef struct {
    SwitchedMatrix  H;
    SwitchedMatrix  Q;
    double         *wr;
    double         *wi;
    int             n;
    int             iters;
    bool            converged;
} QRWorkspace;

/** Switching sequence specification for simulation. */
typedef struct {
    int      *mode_order;
    double   *durations;
    int       length;
    bool      repeat;
} SwitchSequence;

/* ==========================================================================
 * Constructor / Destructor / Utility Declarations
 * ========================================================================== */

/* --- Vector operations (BLAS-like level 1) --- */
SwitchedVector  sv_create(int n);
void            sv_free(SwitchedVector *v);
void            sv_set(SwitchedVector *v, int i, double val);
double          sv_get(const SwitchedVector *v, int i);
double          sv_norm(const SwitchedVector *v);
double          sv_norm_inf(const SwitchedVector *v);
double          sv_dot(const SwitchedVector *a, const SwitchedVector *b);
void            sv_copy(SwitchedVector *dst, const SwitchedVector *src);
void            sv_scale(SwitchedVector *v, double alpha);
void            sv_axpy(SwitchedVector *y, double alpha, const SwitchedVector *x);
void            sv_print(const SwitchedVector *v, FILE *fp);

/* --- Matrix operations --- */
SwitchedMatrix  sm_create(int rows, int cols);
void            sm_free(SwitchedMatrix *M);
void            sm_set(SwitchedMatrix *M, int i, int j, double val);
double          sm_get(const SwitchedMatrix *M, int i, int j);
void            sm_identity(SwitchedMatrix *M);
void            sm_copy(SwitchedMatrix *dst, const SwitchedMatrix *src);
void            sm_transpose(SwitchedMatrix *dst, const SwitchedMatrix *src);
void            sm_add(SwitchedMatrix *C, const SwitchedMatrix *A, const SwitchedMatrix *B);
void            sm_sub(SwitchedMatrix *C, const SwitchedMatrix *A, const SwitchedMatrix *B);
void            sm_mul(SwitchedMatrix *C, const SwitchedMatrix *A, const SwitchedMatrix *B);
void            sm_mul_scalar(SwitchedMatrix *M, double scalar);
void            sm_matvec_mul(SwitchedVector *y, const SwitchedMatrix *A, const SwitchedVector *x);
double          sm_trace(const SwitchedMatrix *M);
double          sm_det_2x2(const SwitchedMatrix *M);
double          sm_det_3x3(const SwitchedMatrix *M);
void            sm_inv_2x2(SwitchedMatrix *inv, const SwitchedMatrix *M);
bool            sm_is_symmetric(const SwitchedMatrix *M);
bool            sm_is_positive_definite(const SwitchedMatrix *M);
double          sm_frobenius_norm(const SwitchedMatrix *M);
void            sm_print(const SwitchedMatrix *M, FILE *fp);

/* --- Lie bracket: [A, B] = AB - BA --- */
void            sm_commutator(SwitchedMatrix *bracket, const SwitchedMatrix *A, const SwitchedMatrix *B);

/* --- Switching signal --- */
SwitchingSignal* ssig_create(int capacity);
void             ssig_free(SwitchingSignal *sig);
void             ssig_record_switch(SwitchingSignal *sig, int new_mode, double time);
void             ssig_set_type(SwitchingSignal *sig, SwitchingSignalType type);
int              ssig_active_mode_at(const SwitchingSignal *sig, double t);
void             ssig_reset(SwitchingSignal *sig);
void             ssig_print(const SwitchingSignal *sig, FILE *fp);

/* --- Subsystem --- */
SwitchedSubsystem* ssub_create(int mode_id, int state_dim);
void               ssub_free(SwitchedSubsystem *sub);
void               ssub_set_hurwitz_check(SwitchedSubsystem *sub);
void               ssub_print(const SwitchedSubsystem *sub, FILE *fp);

/* --- Switched system --- */
SwitchedSystem* sss_create(const char *name, int state_dim, int n_modes);
void            sss_free(SwitchedSystem *sys);
void            sss_add_subsystem(SwitchedSystem *sys, int mode_id, const SwitchedMatrix *A);
void            sss_set_signal(SwitchedSystem *sys, SwitchingSignal *sig);
void            sss_set_initial_state(SwitchedSystem *sys, const SwitchedVector *x0);
void            sss_step(SwitchedSystem *sys, double dt);
void            sss_simulate(SwitchedSystem *sys, const SolverConfig *cfg, const SwitchSequence *seq);
void            sss_update_lyapunov_values(SwitchedSystem *sys);
void            sss_print(const SwitchedSystem *sys, FILE *fp);

/* --- Lie algebra --- */
LieAlgebraCondition* lie_check_create(int n_modes, int dim);
void                 lie_check_free(LieAlgebraCondition *la);
void                 lie_compute_bracket(SwitchedMatrix *bracket, const SwitchedMatrix *A, const SwitchedMatrix *B);
bool                 lie_check_commute(const SwitchedMatrix *A, const SwitchedMatrix *B);
bool                 lie_check_simult_triang(SwitchedMatrix **A, int n_modes, int dim);
bool                 lie_condition_solvable(LieAlgebraCondition *la, SwitchedMatrix **A, int n_modes, int dim);

/* --- Eigenvalue and QR --- */
QRWorkspace* qr_create(int n);
void         qr_free(QRWorkspace *qr);
int          qr_eigenvalues(QRWorkspace *qr, SwitchedMatrix *M, EigenvalueResult *results);
double       spectral_radius(const SwitchedMatrix *M);
bool         is_hurwitz_matrix(const SwitchedMatrix *M);

#endif /* SWITCHED_TYPES_H */