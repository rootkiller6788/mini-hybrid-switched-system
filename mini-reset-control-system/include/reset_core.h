/* reset_core.h - Core Definitions for Reset Control Systems
 *
 * Reset control systems form a subclass of hybrid systems where
 * the controller state undergoes instantaneous reset (jump) when
 * the input signal crosses zero or satisfies a trigger condition.
 * This nonlinear mechanism can overcome fundamental linear limitations
 * such as Bode's integral constraint and waterbed effect.
 *
 * Knowledge coverage:
 *   L1: Reset controller, reset condition, reset map, after-reset state,
 *       reset memory, trigger types
 *   L2: Reset times, reset interval, reset band, reset ratio, dwell time,
 *       Zeno behavior prevention, chattering avoidance
 *   L3: State-space (A,B,C,D,Ar), flow set F, jump set J,
 *       base-linear-plus-reset hybrid automaton structure
 *
 * Primary references:
 *   [BB12] Banos & Barreiro, "Reset Control Systems" (2012), Springer
 *   [Cle58] Clegg, "A nonlinear integrator for servomechanisms" (1958),
 *           Trans. AIEE, vol.77, pp.41-42
 *   [HR75] Horowitz & Rosenbaum, "Non-linear design for cost of feedback
 *          reduction in systems with large parameter uncertainty" (1975),
 *          Int. J. Control, vol.21, pp.977-1001
 *   [ZCH00] Zheng, Chait, Hollot, "Reset control systems: the zero-crossing
 *           resetting law" (2000), Automatica, vol.36, pp.1861-1867
 *   [NZT08] Nesic, Zaccarian, Teel, "Stability properties of reset systems"
 *           (2008), IEEE Trans. Autom. Control, vol.53, pp.542-549
 *   [GST12] Goebel, Sanfelice, Teel, "Hybrid Dynamical Systems" (2012),
 *           Princeton University Press
 */

#ifndef RESET_CORE_H
#define RESET_CORE_H

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L1: Core Definitions - Reset Control Primitive Enum Types
 * ================================================================ */

/** ResetTriggerType enumerates what type of signal event triggers
 *  the instantaneous reset of the controller state.
 *
 *  The trigger condition determines the jump set J in hybrid formalism:
 *    J = { (x,e) : trigger condition holds }
 *
 *  Ref: [BB12] Section 2.2; [GST12] Chapter 1
 */
typedef enum {
    RESET_TRIGGER_ZERO_CROSSING = 0,  /* e(t)*e(t^-) <= 0 (ZC law)   */
    RESET_TRIGGER_RISING        = 1,  /* rising ZC only: e crosses up */
    RESET_TRIGGER_FALLING       = 2,  /* falling ZC only: e crosses down */
    RESET_TRIGGER_PERIODIC      = 3,  /* reset at fixed interval Ts   */
    RESET_TRIGGER_STATE_COND    = 4,  /* reset when g(x_c) <= threshold */
    RESET_TRIGGER_TIME_DEP      = 5   /* reset at specified times {tk} */
} ResetTriggerType;

/** ResetMapType enumerates what transformation is applied to the
 *  controller state at the reset instant.
 *
 *  x_c^+ = f_reset(x_c^-, e, t)
 *
 *  Ref: [BB12] Section 2.3, Table 2.1
 */
typedef enum {
    RESET_MAP_ZERO      = 0,  /* x_c^+ = 0                (full reset) */
    RESET_MAP_RATIO     = 1,  /* x_c^+ = rho * x_c^-      (proportional) */
    RESET_MAP_MATRIX    = 2,  /* x_c^+ = Ar * x_c^-       (matrix) */
    RESET_MAP_SETPOINT  = 3,  /* x_c^+ = x_sp             (prestored) */
    RESET_MAP_FUNCTION  = 4   /* x_c^+ = phi(x_c^-, e, t) (general) */
} ResetMapType;

/** ResetResult enumerates possible outcomes of a reset operation.
 *  Used as return codes from reset_check_and_apply() and similar. */
typedef enum {
    RESET_OK            = 0,  /* reset applied successfully         */
    RESET_NO_CROSSING   = 1,  /* no zero crossing detected          */
    RESET_NOT_ARMED     = 2,  /* crossing signal outside reset band */
    RESET_DWELL_BLOCKED = 3,  /* dwell time not yet elapsed         */
    RESET_CHATTERED     = 4,  /* potential Zeno behavior detected   */
    RESET_NULL_PTR      = 5   /* invalid NULL pointer argument      */
} ResetResult;

/* ================================================================
 * L2: Core Concepts - Reset Dynamics Parameter Structures
 * ================================================================ */

/** ResetBand - tolerance band around zero crossing for avoiding
 *  infinite chattering (Zeno behavior).
 *
 *  The band defines a region [-U, -L] U [L, U] where:
 *    - When |e| > band_upper: reset is disarmed (controller far from ZC)
 *    - When |e| < band_lower: reset is armed (controller near ZC)
 *    - band_hyst provides hysteresis to prevent rapid arm/disarm cycling
 *
 *  Ref: [BB12] Section 4.2.1 on chattering avoidance
 */
typedef struct {
    double band_lower;   /**< lower threshold to arm reset, > 0      */
    double band_upper;   /**< upper threshold to disarm, > band_lower */
    double band_hyst;    /**< hysteresis band to prevent chattering   */
} ResetBand;

/** ResetRatio - the fraction of the controller state that survives
 *  after reset: x_c^+ = rho * x_c^-.
 *
 *  For Clegg integrator: rho = 0 (full reset, complete state erasure).
 *  For FORE: 0 < rho < 1 (partial reset, fractional state preservation).
 *  rho = 1 corresponds to no reset (purely linear behavior).
 *
 *  Time-varying reset ratio rho(t) can provide additional design
 *  freedom for trade-off between disturbance rejection and noise.
 *
 *  Ref: [HR75] on FORE; [BB12] Section 3.3 on variable reset ratio
 */
typedef struct {
    double     rho;              /**< nominal reset ratio in [0, 1)  */
    double     rho_min;          /**< minimum ratio if time-varying  */
    double     rho_max;          /**< maximum ratio if time-varying  */
    bool       is_time_varying;  /**< true if rho varies with time   */
} ResetRatio;

/** ResetInterval - accumulated statistics of reset timing events.
 *
 *  Useful for detecting Zeno phenomena (infinitely many resets in
 *  finite time) and for tuning dwell-time parameters.
 *
 *  Zeno condition: if t_min_interval -> 0 as t increases, the
 *  system exhibits Zeno behavior and the hybrid model is degenerate.
 *
 *  Ref: [BB12] Section 4.3 on reset timing analysis
 */
typedef struct {
    double     t_last;           /**< time of last reset event       */
    double     t_min_interval;   /**< minimum observed reset interval*/
    double     t_max_interval;   /**< maximum observed reset interval*/
    double     t_avg_interval;   /**< average reset interval (EMA)   */
    int        n_resets;         /**< cumulative number of resets    */
} ResetInterval;

/** ResetTimeSequence - predefined schedule of reset instants.
 *
 *  For time-dependent reset trigger type:
 *    t_k in {t_0, t_1, ..., t_{N-1}}
 *
 *  The idx_next pointer advances after each reset to the next
 *  scheduled instant. Useful for periodic or aperiodic reset.
 *
 *  Ref: [ZCH00] on periodic reset control
 */
typedef struct {
    double    *times;            /**< ordered array of reset instants  */
    int        n_times;          /**< number of scheduled resets       */
    int        idx_next;         /**< index of next reset to trigger   */
} ResetTimeSequence;

/* ================================================================
 * L1 continued: Core Structs - System Components
 * ================================================================ */

/** ResetLinearBase - continuous-time linear base system describing
 *  the flow dynamics between reset instants.
 *
 *  Flow equation (on flow set C, when no reset occurs):
 *    dx/dt = A x + B u
 *    y     = C x + D u
 *
 *  where x in R^n (state), u in R^m (input), y in R^p (output).
 *  Matrices stored in row-major order for C compatibility.
 *
 *  For discrete-time systems:
 *    x[k+1] = A x[k] + B u[k]
 *    y[k]   = C x[k] + D u[k]
 *
 *  Ref: [BB12] Equation (2.1); [GST12] Section 2.1
 */
typedef struct {
    int        n;     /**< state dimension                           */
    int        m;     /**< input dimension                           */
    int        p;     /**< output dimension                          */
    double    *A;     /**< system matrix, row-major, size n*n        */
    double    *B;     /**< input matrix, row-major, size n*m         */
    double    *C;     /**< output matrix, row-major, size p*n        */
    double    *D;     /**< feedthrough matrix, row-major, size p*m   */
} ResetLinearBase;

/** ResetJumpMap - the reset (jump) map applied instantaneously
 *  at each reset event.
 *
 *  Jump equation (on jump set D, when reset triggers):
 *    x_c^+ = Ar * x_c^- + Br * e
 *
 *  where Ar in R^{nc*nc} is the reset matrix and Br in R^{nc*1}
 *  optionally feeds the error signal into the after-reset state.
 *
 *  Key special cases (see ResetMapType for semantic classification):
 *    Ar = 0_nc,      Br = 0     ==> full reset (Clegg integrator)
 *    Ar = rho * I,   Br = 0     ==> proportional reset (FORE)
 *    Ar = diag(a_i), Br = 0     ==> partial reset (selective states)
 *    Ar = I,         Br != 0    ==> additive reset (error injection)
 *
 *  Ref: [BB12] Equation (2.5); [NZT08] Section 2
 */
typedef struct {
    double    *Ar;    /**< reset matrix, row-major, nc*nc            */
    double    *Br;    /**< reset input vector, length nc             */
    int        nc;    /**< controller state dimension                 */
} ResetJumpMap;

/** ResetCondition - algebraic condition defining when the reset
 *  (jump) should occur, i.e., the jump set D.
 *
 *  The condition is evaluated at each time step (discrete) or
 *  via event detection (continuous simulation):
 *
 *  For zero-crossing triggers:
 *    Reset if e(t) * e(t^-) <= 0  and  |e(t^-)| > epsilon
 *
 *  For state-condition triggers:
 *    Reset if g(x_c) = sum_i H[i] * x_c[i] <= threshold
 *
 *  t_min_dwell enforces a minimum time between consecutive resets
 *  to guarantee non-Zeno behavior (dwell-time regularization).
 *
 *  Ref: [BB12] Section 4.2; [GST12] Definition 2.10
 */
typedef struct {
    ResetTriggerType  trig_type;    /**< type of trigger condition    */
    double           *H;            /**< condition vector (state trig)*/
    double            threshold;    /**< threshold for state condition*/
    double            epsilon;      /**< numerical zero-cross tolerance*/
    double            t_min_dwell;  /**< minimum dwell between resets */
} ResetCondition;

/* ================================================================
 * L3: Mathematical Structures - Complete Reset System
 * ================================================================ */

/** ResetSystem - Complete reset control system descriptor that
 *  bundles the flow dynamics, jump map, trigger condition, and
 *  runtime state into a single structure.
 *
 *  The hybrid automaton H = (C, F, D, G) is defined by:
 *    - Flow set C = R^n \ D  (or a subset for banded reset)
 *    - Flow map  F(x,u) = A x + B u
 *    - Jump set  D = { x : trigger condition holds }
 *    - Jump map  G(x,e) = Ar x + Br e
 *
 *  The solution concept follows Goebel, Sanfelice, Teel (2012):
 *  a hybrid arc defined on a hybrid time domain (t,j) where
 *  t is the continuous time and j is the number of jumps.
 *
 *  Runtime fields (x_c, t, e_prev, reset_armed) are mutated by
 *  the simulation engine and reset-check functions.
 *
 *  Ref: [GST12] Chapter 2; [BB12] Section 2.4
 */
typedef struct {
    ResetLinearBase    *flow;        /**< continuous flow dynamics    */
    ResetJumpMap       *jump;        /**< discrete jump map           */
    ResetCondition     *cond;        /**< reset trigger condition     */
    ResetRatio         *ratio;       /**< reset ratio parameters      */
    ResetBand          *band;        /**< zero-cross detection band   */
    ResetInterval      *interval;    /**< accumulated reset statistics*/
    int                 nc;          /**< controller state dimension  */
    double             *x_c;         /**< current controller state    */
    double              t;           /**< current simulation time     */
    bool                dt_mode;     /**< true if discrete-time mode  */
    double              ts;          /**< sample time (dt_mode)       */
    bool                reset_armed; /**< armed for next ZC trigger   */
    double              e_prev;      /**< previous error for ZC detect*/
} ResetSystem;

/** ResetMemory - stores the before-and-after reset state snapshots
 *  for post-simulation analysis, stability verification, and
 *  reset-to-reset mapping studies.
 *
 *  Used in the reset-to-reset analysis framework [BB12, Ch.5],
 *  which studies the discrete-time map from one reset event
 *  to the next, forming a Poincare-like analysis for hybrid systems.
 */
typedef struct {
    double     *x_before;     /**< controller state right before reset*/
    double     *x_after;      /**< controller state right after reset */
    double      t_reset;      /**< time instant of this reset        */
    double      e_at_reset;   /**< error signal value at reset moment*/
    int         reset_seq;    /**< sequence number (0, 1, 2, ...)    */
} ResetMemory;

/* ================================================================
 * L2: Core Operations - Lifecycle and Configuration Functions
 * ================================================================ */

/** Allocate and zero-initialize a linear base system of given
 *  dimensions (n states, m inputs, p outputs).
 *  Complexity: O(n^2 + n*m + p*n + p*m) memory, O(1) time.
 *  Returns NULL on dimension error or allocation failure. */
ResetLinearBase* reset_base_create(int n, int m, int p);

/** Deep-copy a linear base system including all matrix data.
 *  Complexity: O(n^2 + n*m + p*n + p*m) time and memory.
 *  Returns NULL if base is NULL or allocation fails. */
ResetLinearBase* reset_base_clone(const ResetLinearBase *base);

/** Free all memory associated with a linear base system.
 *  Safe to call with NULL (no-op). */
void reset_base_free(ResetLinearBase *base);

/** Allocate a complete reset system descriptor with nc controller
 *  states. All sub-structures are allocated and zero-initialized.
 *  The flow, jump, and cond sub-structures are allocated but
 *  their internal matrices must be set separately via setter functions.
 *  Returns NULL on allocation failure or nc <= 0. */
ResetSystem* reset_sys_create(int nc);

/** Free a reset system and all sub-structures recursively.
 *  Sets pointers to NULL after freeing. Safe with NULL input. */
void reset_sys_free(ResetSystem *rsys);

/** Create a standard Clegg integrator [Cle58]:
 *  Flow: dx/dt = e (integrator), y = x
 *  Jump: x^+ = 0 (full reset on zero crossing)
 *  State-space: A=0, B=1, C=1, D=0; Ar=0, Br=0.
 *  Trigger: zero-crossing (any sign change).
 *  Returns fully initialized reset system. */
ResetSystem* reset_clegg_create(void);

/** Create a FORE (First Order Reset Element) [HR75]:
 *  Base transfer function: G(s) = K / (tau*s + 1)
 *  Flow: A = -1/tau, B = K/tau, C = 1, D = 0
 *  Jump: x^+ = rho * x^- (proportional reset)
 *  Trigger: zero-crossing (any sign change).
 *  Parameters: K = DC gain, tau = time constant, rho in [0,1).
 *  Returns NULL if tau <= 0 or rho < 0 or rho >= 1. */
ResetSystem* reset_fore_create(double K, double tau, double rho);

/** Set the reset ratio for a proportional-reset element.
 *  Updates rsys->ratio->rho and syncs Ar = rho * I.
 *  0 <= rho < 1. rho=0 is full reset.
 *  Returns RESET_NULL_PTR if rsys is NULL, RESET_OK otherwise.
 *  Complexity: O(nc^2) to update the Ar matrix. */
ResetResult reset_set_ratio(ResetSystem *rsys, double rho);

/** Evaluate the reset trigger condition and apply the jump map
 *  if the condition is satisfied.
 *
 *  For zero-crossing triggers:
 *    - Checks if sign(e) != sign(e_prev) && |e_prev| > cond->epsilon
 *    - Respects reset band arming/disarming logic
 *    - Enforces dwell time constraint
 *    - Updates x_c, e_prev, t, and interval statistics
 *
 *  For time-dependent triggers:
 *    - Checks if current time t has crossed next scheduled reset time
 *    - Advances idx_next and applies jump map
 *
 *  Complexity: O(nc^2) for matrix-vector multiply at reset.
 *  Returns RESET_OK if reset was applied, error code otherwise. */
ResetResult reset_check_and_apply(ResetSystem *rsys, double e, double e_prev);

/** Apply a manual (forced) reset at given time, bypassing all
 *  trigger conditions, band checks, and dwell time constraints.
 *  Stores before/after state in a ResetMemory snapshot.
 *  Complexity: O(nc^2). Returns RESET_NULL_PTR if rsys is NULL. */
ResetResult reset_apply_manual(ResetSystem *rsys, double t, double e);

/** Configure zero-crossing trigger with reset band parameters.
 *  trig_type: which type of ZC to detect (any, rising, falling)
 *  band_lower: arming threshold (must be > 0)
 *  band_upper: disarming threshold (must be > band_lower)
 *  hyst: hysteresis for arm/disarm transition
 *  Complexity: O(1). */
void reset_config_zc_trigger(ResetSystem *rsys, ResetTriggerType trig_type,
                              double band_lower, double band_upper,
                              double hyst);

/** Set the minimum dwell time between consecutive resets.
 *  This is the primary mechanism for preventing Zeno behavior.
 *  t_min_dwell > 0 guarantees that the number of resets in any
 *  finite time interval is bounded.
 *  Ref: [NZT08] Theorem 1 on dwell-time regularization. */
void reset_set_dwell_time(ResetSystem *rsys, double t_min_dwell);

/** Configure a time-dependent reset schedule.
 *  The times array must be sorted in strictly increasing order.
 *  Makes an internal copy of the times array.
 *  Complexity: O(n) time and memory.
 *  Ref: [ZCH00] on periodic reset control design. */
void reset_set_time_schedule(ResetSystem *rsys, const double *times, int n);

/** Retrieve the accumulated reset interval statistics (read-only).
 *  Returns pointer to internal ResetInterval struct.
 *  Returns NULL if rsys is NULL. */
const ResetInterval* reset_get_interval_stats(const ResetSystem *rsys);

/** Set the linear base system matrices from raw arrays.
 *  Takes ownership of the pointer arrays (they will be freed
 *  by reset_sys_free). The caller should not free them afterwards.
 *  Old matrices in rsys->flow are freed automatically.
 *  Complexity: O(1). */
void reset_set_linear_base(ResetSystem *rsys, int n, int m, int p,
                            double *A, double *B, double *C, double *D);

/** Set the jump map matrices from raw arrays.
 *  Takes ownership of the pointer arrays.
 *  Complexity: O(1). */
void reset_set_jump_map(ResetSystem *rsys, double *Ar, double *Br, int nc);

/** Evaluate the flow derivative dx/dt = A x + B u.
 *  dxdt must be pre-allocated with n elements.
 *  Complexity: O(n^2 + n*m).
 *  Returns 0 on success, -1 on NULL pointer. */
int reset_eval_flow_deriv(const ResetLinearBase *base,
                           const double *x, const double *u, double *dxdt);

/** Evaluate the output equation y = C x + D u.
 *  y must be pre-allocated with p elements.
 *  Complexity: O(p*n + p*m).
 *  Returns 0 on success, -1 on NULL pointer. */
int reset_eval_output(const ResetLinearBase *base,
                       const double *x, const double *u, double *y);

#ifdef __cplusplus
}
#endif

#endif /* RESET_CORE_H */