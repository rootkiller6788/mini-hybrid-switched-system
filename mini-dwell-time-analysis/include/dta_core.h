#ifndef DTA_CORE_H
#define DTA_CORE_H
#include <stdbool.h>
#include <stddef.h>

/* ==============================================================
 * dta_core.h - Dwell-Time Analysis Core Types
 *
 * Dwell-time analysis studies stability of switched systems
 * where the switching signal σ(t) has a minimum dwell time τ_d
 * between consecutive switches.
 *
 * A switched system is given by:
 *   ẋ(t) = f_{σ(t)}(x(t)),  σ(t) ∈ {1, 2, ..., m}
 *
 * where σ(t) is a piecewise-constant switching signal.
 *
 * Key quantities:
 *   τ_d = inf_{k} (t_{k+1} - t_k)  — dwell time
 *   τ_a = average dwell time (over long intervals)
 *   N_σ(T,t) = number of switches in [t, T]
 *
 * Fundamental result (Morse 1996, Hespanha & Morse 1999):
 *   If all subsystems are stable and τ_d is sufficiently large,
 *   the switched system is globally asymptotically stable (GAS).
 *
 * References:
 *   Liberzon (2003) "Switching in Systems and Control"
 *   Liberzon & Morse (1999) "Basic problems in stability and
 *     design of switched systems", IEEE Control Systems
 *   Morse (1996) "Supervisory control of families of linear
 *     set-point controllers", IEEE TAC
 *   Hespanha & Morse (1999) "Stability of switched systems
 *     with average dwell-time", IEEE CDC
 * ============================================================== */

/* --- Mode types --- */
typedef enum {
    DTA_MODE_STABLE = 0,
    DTA_MODE_UNSTABLE = 1,
    DTA_MODE_MARGINAL = 2,
    DTA_MODE_UNKNOWN = 3
} DTA_ModeStability;

/* --- System class --- */
typedef enum {
    DTA_LINEAR = 0,
    DTA_AFFINE = 1,
    DTA_NONLINEAR = 2,
    DTA_BILINEAR = 3
} DTA_SystemClass;

/* --- Stability result --- */
typedef enum {
    DTA_GAS = 0,
    DTA_GES = 1,
    DTA_GUAS = 2,
    DTA_GUES = 3,
    DTA_UNSTABLE = 4,
    DTA_INCONCLUSIVE = 5,
    DTA_STABLE_NOT_ASYM = 6
} DTA_StabilityVerdict;

/* --- Dwell time characterization --- */
typedef enum {
    DTA_DWELL_CONSTANT = 0,
    DTA_DWELL_AVERAGE = 1,
    DTA_DWELL_MODE_DEP = 2,
    DTA_DWELL_PERSISTENT = 3,
    DTA_DWELL_RANGE = 4
} DTA_DwellType;

/**
 * Single mode of a switched system.
 * For linear modes:  ẋ = A x + B u
 * For nonlinear modes: ẋ = f(x, u) via function pointer
 */
typedef struct {
    int id;
    double* A;
    double* B;
    double* C;
    int n;
    int p;
    int q;
    bool has_affine_term;
    double* affine_b;
    DTA_ModeStability stability;
    double max_eigenvalue_real;
    double spectral_radius;
    void (*dynamics)(double t, const double* x, int n, const double* u,
                     int p, double* dx, void* context);
    void* context;
} DTA_SystemMode;

/**
 * A switched system S = {f_1, f_2, ..., f_m} with m modes.
 */
typedef struct {
    DTA_SystemMode* modes;
    int n_modes;
    int state_dim;
    int input_dim;
    int output_dim;
    DTA_SystemClass class;
    double current_time;
    int current_mode;
    int max_switches;
} DTA_SwitchedSystem;

/**
 * Dwell time configuration.
 */
typedef struct {
    DTA_DwellType type;
    double tau_d;
    double tau_d_min;
    double tau_d_max;
    double* mode_tau;
    double N0;
    double T;
    double t0;
} DTA_DwellTimeConfig;

/**
 * A switching signal σ: [t0, T] → {1, ..., m}.
 */
typedef struct {
    double* switch_times;
    int* mode_sequence;
    int n_switches;
    int capacity;
    double t_start;
    double t_end;
} DTA_SwitchingSignal;

/**
 * State trajectory of a switched system.
 */
typedef struct {
    double* t_samples;
    double** x_samples;
    int n_samples;
    int capacity;
    int state_dim;
} DTA_StateTrajectory;

/**
 * Solution of the Lyapunov equation A^T P + P A = -Q for a mode.
 */
typedef struct {
    double* P;
    double* Q;
    int n;
    int mode_id;
    double min_eig_P;
    double max_eig_P;
    double min_eig_Q;
    bool solved;
} DTA_LyapunovSolution;

/* --- Core API --- */

DTA_SwitchedSystem* dta_system_create(int n_modes, int state_dim,
                                       int input_dim, int output_dim);
void dta_system_free(DTA_SwitchedSystem* sys);
int dta_system_set_linear_mode(DTA_SwitchedSystem* sys, int mode_idx,
                                const double* A, const double* B,
                                const double* C);
int dta_system_set_nonlinear_mode(DTA_SwitchedSystem* sys, int mode_idx,
    void (*dynamics)(double, const double*, int, const double*, int,
                     double*, void*), void* context);
int dta_system_set_affine_mode(DTA_SwitchedSystem* sys, int mode_idx,
                                const double* A, const double* b);
void dta_system_rhs(const DTA_SwitchedSystem* sys, int mode,
                    double t, const double* x, const double* u,
                    double* dx);
int dta_eigenvalues(const double* A, int n, double* real_parts,
                    double* imag_parts, int max_iter);
bool dta_is_hurwitz(const double* A, int n, double tol);
bool dta_is_schur(const double* A, int n, double tol);
void dta_matrix_exp(const double* A, int n, double t, double* result);
double dta_matrix_measure(const double* A, int n);
int dta_solve_lyapunov(const double* A, int n, const double* Q,
                        double* P);
double dta_l2_gain(const double* A, const double* B, const double* C,
                    int n, int p, int q);

#endif /* DTA_CORE_H */
