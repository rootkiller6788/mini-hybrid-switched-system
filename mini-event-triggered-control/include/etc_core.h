#ifndef ETC_CORE_H
#define ETC_CORE_H

#include <stdbool.h>
#include <stddef.h>

/* ============================================================================
 * Event-Triggered Control (ETC) — Core Definitions
 *
 * Event-triggered control replaces periodic sampling with state-dependent
 * triggering. Control updates occur only when a triggering condition is
 * violated, reducing communication/computation while preserving stability.
 *
 * Foundational references:
 *   Tabuada (2007) — IEEE TAC: Event-triggered real-time scheduling
 *   Heemels, Johansson, Tabuada (2012) — An introduction to ETC & STC
 *   Astrom & Bernhardsson (2002) — Comparison of Riemann & Lebesgue sampling
 *   Lemmon (2010) — Event-triggered feedback in control, estimation, optimization
 * ============================================================================ */

/* --- Fundamental Enumerations --- */

/** Trigger mechanism classification.
 *  STATIC: threshold depends only on current state/error
 *  DYNAMIC: threshold evolves with an internal dynamic variable
 *  MIXED: combines state-dependent and time-dependent conditions
 *  ABSOLUTE: fixed threshold (not state-relative)
 *  SELF_TRIGGERED: next event time computed at current event */
typedef enum {
    ETC_TRIGGER_STATIC = 0,
    ETC_TRIGGER_DYNAMIC = 1,
    ETC_TRIGGER_MIXED = 2,
    ETC_TRIGGER_ABSOLUTE = 3,
    ETC_TRIGGER_SELF_TRIGGERED = 4
} ETCTriggerType;

/** System regime classification during event-triggered operation.
 *  ASYMPTOTIC_STABLE: system converges to origin
 *  ISS_PRACTICAL: input-to-state practically stable
 *  ZENO: infinite events in finite time (pathological)
 *  MINIMAL_INTERVAL: provably positive minimum inter-event time
 *  UNKNOWN: regime not yet classified */
typedef enum {
    ETC_REGIME_ASYMPTOTIC = 0,
    ETC_REGIME_ISS_PRACTICAL = 1,
    ETC_REGIME_ZENO = 2,
    ETC_REGIME_MINIMAL_INTERVAL = 3,
    ETC_REGIME_UNKNOWN = 4
} ETCRegime;

/** Comparison function classes (K, K∞, KL functions).
 *  Used in ISS Lyapunov characterizations. */
typedef enum {
    ETC_K_CLASS = 0,      /* α: R≥0 → R≥0, α(0)=0, strictly increasing */
    ETC_KINF_CLASS = 1,   /* K-class + lim_{s→∞} α(s) = ∞ */
    ETC_KL_CLASS = 2      /* β: R≥0 × R≥0 → R≥0, β(·,t) ∈ K, β(s,·) decreasing → 0 */
} ETCComparisonClass;

/** Sampling mechanism type */
typedef enum {
    ETC_SAMPLING_RIEMANN = 0,  /* Periodic sampling (time-triggered) */
    ETC_SAMPLING_LEBESGUE = 1  /* Event-triggered sampling (state-triggered) */
} ETCSamplingType;

/* --- Core Numeric Types --- */

/** N-dimensional vector (state, error, or control). */
typedef struct {
    double* data;      /* Component array */
    int dim;           /* Dimension */
} ETCVector;

/** Dense matrix for state-space representation. */
typedef struct {
    double* data;      /* Row-major flat array */
    int rows;
    int cols;
} ETCMatrix;

/** Quadratic Lyapunov function V(x) = x^T P x.
 *  P ∈ R^{n×n}, symmetric positive definite. */
typedef struct {
    ETCMatrix P;            /* Lyapunov matrix */
    double P_norm;          /* ||P|| for quick reference */
    double lambda_min_P;    /* Minimum eigenvalue of P */
    double lambda_max_P;    /* Maximum eigenvalue of P */
    bool is_positive_definite;
} ETCLyapunovFunction;

/* --- Event Record --- */

/** A single event in the control timeline.
 *  Records the state at which control was updated and the triggering condition
 *  value at that instant. */
typedef struct {
    double time;               /* Event occurrence time */
    int event_index;           /* Sequential event number (0 = initial) */
    ETCVector state;           /* Sampled state x(t_k) */
    ETCVector control;         /* Computed control u(t_k) */
    ETCVector error_at_event;  /* Measurement error e(t_k) = 0 at event */
    double trigger_value;      /* Γ(x,e) at event instant */
    double inter_event_time;   /* Δ_k = t_k - t_{k-1} */
} ETCEvent;

/** Event history buffer for post-simulation analysis. */
typedef struct {
    ETCEvent* events;
    int n_events;
    int capacity;
    double min_iet;          /* Minimum inter-event time observed */
    double max_iet;          /* Maximum inter-event time observed */
    double avg_iet;          /* Average inter-event time */
    bool zeno_detected;      /* True if infinite events in finite window */
} ETCEventHistory;

/* --- Core ETC System Model --- */

/** Complete event-triggered control system.
 *
 *  Plant:  ẋ(t) = A x(t) + B u(t)   (continuous-time linear system)
 *  Controller: u(t) = K x(t_k),  t ∈ [t_k, t_{k+1})
 *  Measurement error: e(t) = x(t_k) − x(t)
 *  Trigger condition: t_{k+1} = inf{ t > t_k : Γ(x(t), e(t)) ≥ 0 }
 *
 *  Closed loop: ẋ(t) = (A + BK) x(t) − BK e(t)
 */
typedef struct {
    /* System matrices */
    ETCMatrix A;             /* State matrix (n×n) */
    ETCMatrix B;             /* Input matrix (n×m) */
    ETCMatrix K;             /* State-feedback gain (m×n) */
    ETCMatrix Acl;           /* Closed-loop: A + BK (n×n), precomputed */
    int n_states;            /* State dimension n */
    int n_inputs;            /* Input dimension m */
    bool is_stabilizable;    /* (A,B) stabilizable check result */

    /* Current state and control */
    ETCVector x;             /* Current continuous state x(t) */
    ETCVector x_hat;         /* Last sampled state x(t_k) */
    ETCVector e;             /* Measurement error e(t) = x_hat − x */
    ETCVector u;             /* Current control input u(t) */

    /* Trigger configuration */
    ETCTriggerType trigger_type;
    double sigma;            /* Relative threshold parameter σ */
    double epsilon;          /* Absolute threshold parameter ε */
    double (*trigger_fn)(const ETCVector* x, const ETCVector* e,
                         double sigma, double epsilon);  /* Γ(x,e) */

    /* Dynamic trigger state (for DYNAMIC type) */
    double eta;              /* Internal dynamic variable η(t) */
    double eta_dot;          /* Derivative of η */
    double theta;            /* Dynamic threshold parameter θ */

    /* Lyapunov function for stability analysis */
    ETCLyapunovFunction V;

    /* Timing */
    double t;                /* Current simulation time */
    double t_last_event;     /* Time of last event */
    double dt;               /* Integration step size */
    double t_max;            /* Simulation horizon */
    int event_count;         /* Total number of events triggered */
    double min_iet_observed; /* Running minimum inter-event time */

    /* Regime classification */
    ETCRegime regime;
    double zeno_bound;       /* Estimated Zeno time (if ZENO detected) */

    /* History */
    ETCEventHistory history;
} ETCSystem;

/* --- Core API --- */

/** Create an ETC system from given matrices and gain.
 *  Allocates and initializes all internal structures.
 *  The trigger function defaults to quadratic state-relative:
 *  Γ(x,e) = |e|² − σ|x|²
 *
 *  @param A  State matrix (n×n), row-major
 *  @param B  Input matrix (n×m), row-major
 *  @param K  Feedback gain (m×n), row-major
 *  @param n  State dimension
 *  @param m  Input dimension
 *  @return   Newly allocated ETCSystem (must be freed with etc_system_free) */
ETCSystem* etc_system_create(const double* A, const double* B,
                              const double* K, int n, int m);

/** Free all memory associated with an ETC system. */
void etc_system_free(ETCSystem* sys);

/** Set the initial state x(0) and reset sampled state x̂(0) = x(0). */
void etc_system_set_initial_state(ETCSystem* sys, const double* x0);

/** Set the trigger function. If NULL, uses default quadratic trigger.
 *  Trigger signature: Γ(x, e, σ, ε) >= 0 ⇒ fire event. */
void etc_system_set_trigger(ETCSystem* sys, ETCTriggerType type,
                             double sigma, double epsilon,
                             double (*trigger_fn)(const ETCVector*, const ETCVector*,
                                                  double, double));

/** Set dynamic trigger parameters (η₀, θ). Only valid for DYNAMIC type. */
void etc_system_set_dynamic_trigger(ETCSystem* sys, double eta0, double theta);

/** Retrieve the current measurement error norm. */
double etc_system_error_norm(const ETCSystem* sys);

/** Retrieve the current state norm. */
double etc_system_state_norm(const ETCSystem* sys);

/** Evaluate the trigger function at current state. Returns Γ(x,e). */
double etc_system_eval_trigger(const ETCSystem* sys);

/** Check if a new event should be triggered: Γ(x,e) ≥ 0. */
bool etc_system_should_trigger(const ETCSystem* sys);

/** Reset the sampled state after an event: x̂ ← x, e ← 0, record event. */
void etc_system_trigger_event(ETCSystem* sys);

/** Compute the control input u = K x̂. Updates internal u vector. */
void etc_system_compute_control(ETCSystem* sys);

/** Compute the state derivative ẋ = A x + B u. Writes to dx (n-elements). */
void etc_system_compute_derivative(const ETCSystem* sys, double* dx);

/** Advance the continuous-time state by one RK4 integration step.
 *  Returns false if the step could not be completed (e.g., Zeno detected). */
bool etc_system_step_rk4(ETCSystem* sys);

/** Advance the dynamic trigger variable η by one Euler step.
 *  η̇ = −β η + σ|x|² − |e|²  (β > 0 ensures non-negative η). */
void etc_system_step_dynamic_trigger(ETCSystem* sys, double beta);

/** Run full simulation from current time to t_max.
 *  Integrates plant dynamics, evaluates trigger at each step,
 *  fires events and records history. */
void etc_system_simulate(ETCSystem* sys, double t_max, double dt);

/** Print a summary of the ETC system state and event history. */
void etc_system_print(const ETCSystem* sys);

/* --- Vector and Matrix Utilities --- */

/** Create a zero-initialized vector of given dimension. */
ETCVector etc_vector_create(int dim);

/** Free a vector's internal data. */
void etc_vector_free(ETCVector* v);

/** Compute Euclidean norm ||v||₂. */
double etc_vector_norm(const ETCVector* v);

/** Compute dot product v1 · v2. Dimensions must match. */
double etc_vector_dot(const ETCVector* v1, const ETCVector* v2);

/** v_out = v1 − v2 (elementwise subtraction). */
void etc_vector_sub(const ETCVector* v1, const ETCVector* v2, ETCVector* v_out);

/** v_out = s * v (scalar multiply). */
void etc_vector_scale(double s, const ETCVector* v, ETCVector* v_out);

/** v_out = v1 + v2 (elementwise addition). */
void etc_vector_add(const ETCVector* v1, const ETCVector* v2, ETCVector* v_out);

/** Copy src vector to dst. Both must have same dimension. */
void etc_vector_copy(const ETCVector* src, ETCVector* dst);

/** Create a zero matrix of given dimensions. */
ETCMatrix etc_matrix_create(int rows, int cols);

/** Free a matrix's internal data. */
void etc_matrix_free(ETCMatrix* m);

/** Matrix-vector multiply: y = A x. Dimensions must be compatible. */
void etc_matrix_vec_mul(const ETCMatrix* A, const ETCVector* x, ETCVector* y);

/** Matrix multiply: C = A × B. Dimensions: A(m×k) × B(k×n) → C(m×n). */
void etc_matrix_mul(const ETCMatrix* A, const ETCMatrix* B, ETCMatrix* C);

/** Matrix add: C = A + B. All same dimensions. */
void etc_matrix_add(const ETCMatrix* A, const ETCMatrix* B, ETCMatrix* C);

/** Scale matrix: B = s * A. */
void etc_matrix_scale(double s, const ETCMatrix* A, ETCMatrix* B);

/** Compute eigenvalues of a 2×2 real matrix.
 *  Returns (λ₁, λ₂) as complex pair in real/im form.
 *  re[0]=Re(λ₁), im[0]=Im(λ₁), re[1]=Re(λ₂), im[1]=Im(λ₂). */
void etc_matrix_eigenvalues_2x2(const ETCMatrix* A,
                                 double* re, double* im);

/** Compute the spectral radius ρ(A) = max|λᵢ| for a square matrix.
 *  Uses power iteration for n > 2, direct formula for n ≤ 2. */
double etc_matrix_spectral_radius(const ETCMatrix* A);

/** Check if matrix is symmetric. */
bool etc_matrix_is_symmetric(const ETCMatrix* A);

/** Check if a symmetric matrix is positive definite (Cholesky attempt).
 *  Returns true if all leading principal minors positive. */
bool etc_matrix_is_positive_definite(const ETCMatrix* A);

/** Compute the induced 2-norm ||A||₂ = √ρ(AᵀA). */
double etc_matrix_norm_2(const ETCMatrix* A);

/* --- Event History API --- */

/** Initialize event history with given capacity. */
void etc_history_init(ETCEventHistory* hist, int capacity);

/** Append an event to history. Expands capacity if needed. */
void etc_history_append(ETCEventHistory* hist, const ETCEvent* event);

/** Compute statistics (min, max, avg inter-event time) from history. */
void etc_history_compute_stats(ETCEventHistory* hist);

/** Free all memory in event history. */
void etc_history_free(ETCEventHistory* hist);

/** Check for Zeno behavior: infinite events in finite time.
 *  Returns true if event frequency exceeds threshold. */
bool etc_history_check_zeno(const ETCEventHistory* hist, double threshold);

#endif /* ETC_CORE_H */
