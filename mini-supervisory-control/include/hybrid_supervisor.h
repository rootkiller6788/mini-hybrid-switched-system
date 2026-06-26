/**
 * @file hybrid_supervisor.h
 * @brief Hybrid Switched System Supervisor
 *
 * Extends DES supervisory control to hybrid systems where:
 * - Continuous dynamics evolve within each discrete mode
 * - Mode switches are discrete events
 * - The supervisor guards mode transitions to enforce safety
 *
 * This bridges DES theory (Ramadge-Wonham) with hybrid automata,
 * enabling switched system control with formal safety guarantees.
 *
 * References:
 * - Koutsoukos et al. (2000) "Supervisory control of hybrid systems"
 * - Lygeros, Johansson et al. (2012) "Hybrid systems modeling & control"
 * - Tabuada (2009) "Verification and control of hybrid systems"
 *
 * Course mapping:
 * - Berkeley EE222: Nonlinear Systems -> hybrid extensions
 * - MIT 6.832: Underactuated Robotics -> contact switching
 * - ETH 227-0220: Model Reduction -> switched system abstraction
 */

#ifndef HYBRID_SUPERVISOR_H
#define HYBRID_SUPERVISOR_H

#include "des_automaton.h"
#include "supervisor.h"
#include <stddef.h>
#include <stdint.h>

/* ---------------------------------------------------------------------------
 * L1: Hybrid automaton definitions
 * ------------------------------------------------------------------------- */

#define HYB_MAX_MODES        32
#define HYB_MAX_GUARDS       64
#define HYB_MAX_INVARIANTS   64
#define HYB_DIM              6

/** Continuous state vector x in R^n */
typedef struct {
    double  data[HYB_DIM];
    uint8_t dim;
} hyb_state_t;

/**
 * @struct hyb_ode_t
 * @brief Continuous dynamics in a mode: dx/dt = f(x, u)
 *
 * Currently implements linear dynamics: dx/dt = A*x + B*u
 * Can be extended to nonlinear flows.
 */
typedef struct {
    double  A[HYB_DIM][HYB_DIM];
    double  B[HYB_DIM][HYB_DIM];
    double  c[HYB_DIM];        /**< Constant offset (affine term) */
    uint8_t dim;
    uint8_t input_dim;
} hyb_ode_t;

/**
 * @struct hyb_invariant_t
 * @brief Invariant condition: must hold while in a mode.
 *
 * g(x) <= 0  for a set of linear inequalities.
 * If violated, a transition MUST occur (forced transition).
 */
typedef struct {
    double  A[HYB_DIM];   /**< Normal vectors of half-spaces */
    double  b;            /**< Right-hand side */
    uint8_t dim;
} hyb_invariant_t;

/**
 * @struct hyb_guard_t
 * @brief Guard condition: must hold for a switch to be enabled.
 *
 * Similar to invariants but used for enabling transitions.
 */
typedef struct {
    double  A[HYB_DIM];
    double  b;
    uint8_t dim;
} hyb_guard_t;

/**
 * @enum hyb_switch_type_t
 * @brief Types of mode switches.
 */
typedef enum {
    HYB_SWITCH_CONTROLLED = 0,    /**< Supervisor decides */
    HYB_SWITCH_AUTONOMOUS = 1,    /**< Occurs when invariant violated */
    HYB_SWITCH_STOCHASTIC  = 2,   /**< Random with given rate */
} hyb_switch_type_t;

/**
 * @struct hyb_transition_t
 * @brief A transition between hybrid modes.
 */
typedef struct {
    uint16_t    src_mode;
    uint16_t    dst_mode;
    uint16_t    event;           /**< DES event triggering this switch */
    hyb_switch_type_t sw_type;
    hyb_guard_t guard;
    /** Reset map: x' = R*x + r (for state jump on transition) */
    double      reset_R[HYB_DIM][HYB_DIM];
    double      reset_r[HYB_DIM];
    uint8_t     reset_dim;
} hyb_transition_t;

/**
 * @struct hyb_mode_t
 * @brief A single mode of the hybrid automaton.
 */
typedef struct {
    uint16_t        id;
    const char     *name;
    hyb_ode_t       dynamics;
    hyb_invariant_t invariant;
    uint16_t        n_invariants;
    int             is_safe;       /**< 1 if this mode satisfies safety spec */
    int             is_target;     /**< 1 if this mode is a target/goal */
} hyb_mode_t;

/**
 * @struct hybrid_automaton_t
 * @brief Hybrid automaton = DES automaton + continuous dynamics per mode.
 */
typedef struct {
    /** DES layer: discrete event dynamics */
    des_automaton_t  des;

    /** Continuous layer: dynamics per mode */
    hyb_mode_t       modes[HYB_MAX_MODES];
    uint16_t         nmodes;

    /** Transitions between modes */
    hyb_transition_t transitions[HYB_MAX_GUARDS];
    uint16_t         ntrans;

    /** Current continuous state and mode */
    hyb_state_t      x;
    uint16_t         current_mode;

    /** Time */
    double           t;

    /** Name */
    char             name[64];
} hybrid_automaton_t;

/* ---------------------------------------------------------------------------
 * L2: Core hybrid operations
 * ------------------------------------------------------------------------- */

void hybrid_automaton_init(hybrid_automaton_t *H, const char *name);

uint16_t hybrid_add_mode(hybrid_automaton_t *H, const char *mode_name,
                          const hyb_ode_t *dynamics);

void hybrid_set_invariant(hybrid_automaton_t *H, uint16_t mode,
                           const hyb_invariant_t *inv);

int hybrid_add_transition(hybrid_automaton_t *H, uint16_t src, uint16_t dst,
                           uint16_t event, hyb_switch_type_t sw_type,
                           const hyb_guard_t *guard,
                           const double reset_R[HYB_DIM][HYB_DIM],
                           const double *reset_r, uint8_t reset_dim);

int hybrid_guard_satisfied(const hybrid_automaton_t *H, uint16_t trans_idx);

int hybrid_invariant_satisfied(const hybrid_automaton_t *H, uint16_t mode);

void hybrid_step_continuous(hybrid_automaton_t *H, double dt,
                             const double *u, uint8_t u_dim);

int hybrid_execute_event(hybrid_automaton_t *H, uint16_t ev);

void hybrid_reset_state(hybrid_automaton_t *H, uint16_t trans_idx);

uint16_t hybrid_enabled_transitions(const hybrid_automaton_t *H,
                                     uint16_t *enabled, uint16_t cap);

/* ---------------------------------------------------------------------------
 * L5: Hybrid supervisor algorithms
 * ------------------------------------------------------------------------- */

int hybrid_supervisor_safe(const hybrid_automaton_t *H,
                            const supervisor_t *S);

void hybrid_supervisor_synthesize(const hybrid_automaton_t *H,
                                   const des_automaton_t *safety_spec,
                                   supervisor_t *S);

void hybrid_abstraction(const hybrid_automaton_t *H,
                         des_automaton_t *abs_G);

int hybrid_simulate(hybrid_automaton_t *H, supervisor_t *S,
                     double T_end, double dt,
                     hyb_state_t *trajectory, uint16_t max_steps,
                     uint16_t *n_steps);

int hybrid_reachability_check(const hybrid_automaton_t *H,
                               const hyb_state_t *init_set,
                               const hyb_state_t *unsafe_set);

void hybrid_print_mode(const hybrid_automaton_t *H, uint16_t mode, FILE *fp);

void hybrid_print_trajectory(const hyb_state_t *traj, uint16_t n_steps,
                              FILE *fp);

void hybrid_automaton_destroy(hybrid_automaton_t *H);

#endif /* HYBRID_SUPERVISOR_H */
