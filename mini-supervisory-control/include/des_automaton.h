/**
 * @file des_automaton.h
 * @brief Discrete Event System (DES) Automaton — Core Data Structures
 *
 * Defines the foundational automaton model used in Ramadge-Wonham
 * Supervisory Control Theory. An automaton is a 5-tuple:
 *   G = (Q, Σ, δ, q0, Qm)
 * where Q = states, Σ = events, δ = transition function,
 * q0 = initial state, Qm = marked states.
 *
 * References:
 * - Ramadge & Wonham (1989) "The control of discrete event systems"
 * - Cassandras & Lafortune (2008) "Introduction to Discrete Event Systems"
 * - Wonham & Cai (2019) "Supervisory Control of Discrete-Event Systems"
 *
 * Course mapping:
 * - Cambridge 4F3: Nonlinear & Predictive Control
 * - CMU 24-654: Systems Thinking
 * - MIT 6.241J: Dynamic Systems & Control
 * - ETH 227-0216: System Identification
 */

#ifndef DES_AUTOMATON_H
#define DES_AUTOMATON_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

/* ---------------------------------------------------------------------------
 * L1: Core Definitions — event types, state, transition, automaton
 * ------------------------------------------------------------------------- */

/** Maximum number of states and events for a single automaton.
 *  Tunable for large-scale systems; used to enable static allocation. */
#define DES_MAX_STATES     256
#define DES_MAX_EVENTS     128
#define DES_MAX_TRANS      1024

/** Special state value for "undefined" — used in partial transition function. */
#define DES_UNDEF_STATE    ((uint16_t)(-1))

/** Marker for disabled event in supervisor control pattern. */
#define DES_EVENT_DISABLED  0xFF

/**
 * @enum des_event_type_t
 * @brief Classification of events in supervisory control.
 *
 * Σ = Σc ∪ Σu  (controllable vs uncontrollable)
 * Σ = Σo ∪ Σuo (observable vs unobservable)
 *
 * Controllable (Σc): Supervisor can disable these events.
 * Uncontrollable (Σu): Supervisor MUST allow these events (nature, faults).
 * Observable (Σo): Supervisor can see these events occur.
 * Unobservable (Σuo): Supervisor cannot directly observe.
 */
typedef enum {
    DES_EVENT_CONTROLLABLE   = 0x01,
    DES_EVENT_UNCONTROLLABLE = 0x02,
    DES_EVENT_OBSERVABLE     = 0x04,
    DES_EVENT_UNOBSERVABLE   = 0x08,
    DES_EVENT_FORCIBLE       = 0x10  /**< Can be forced to occur (optional) */
} des_event_type_t;

/**
 * @struct des_event_t
 * @brief A single event in the DES alphabet Σ.
 */
typedef struct {
    uint16_t        idx;
    const char     *label;
    uint8_t         flags;   /**< Bitmask of des_event_type_t */
    double          cost;    /**< Optional cost/priority for optimal supervision */
} des_event_t;

/**
 * @struct des_transition_t
 * @brief A single transition in the automaton.
 *
 * Transition: from state `src`, on event `ev`, go to state `dst`.
 * Corresponds to δ: Q × Σ → Q.
 */
typedef struct {
    uint16_t    src;
    uint16_t    ev;
    uint16_t    dst;
} des_transition_t;

/**
 * @struct des_automaton_t
 * @brief Deterministic finite automaton representing a DES plant.
 *
 * G = (Q, Σ, δ, q0, Qm)
 * - Q: finite set of states, numbered 0..nstates-1
 * - Σ: finite alphabet of events
 * - δ: partial transition function Q × Σ → Q
 * - q0: initial state
 * - Qm ⊆ Q: set of marked states (completed tasks)
 *
 * L(G)  = {s ∈ Σ* | δ(q0, s) is defined}  — closed behavior (all strings)
 * Lm(G) = {s ∈ L(G) | δ(q0, s) ∈ Qm}     — marked behavior (completed)
 */
typedef struct {
    /* Σ: alphabet */
    des_event_t  events[DES_MAX_EVENTS];
    uint16_t     nevents;
    uint16_t     n_controllable;
    uint16_t     n_observable;

    /* Q: state set */
    uint16_t     nstates;
    uint16_t     q0;
    uint16_t     n_marked;
    uint16_t     marked[DES_MAX_STATES];

    /* δ: transition function */
    des_transition_t transitions[DES_MAX_TRANS];
    uint16_t     ntrans;

    /** δ adjacency matrix: delta[src][ev] = dst or DES_UNDEF_STATE */
    uint16_t     delta[DES_MAX_STATES][DES_MAX_EVENTS];

    /** Name for debugging */
    char         name[64];
} des_automaton_t;

/* ---------------------------------------------------------------------------
 * L2: Core Concepts — State accessors, event classification, language ops
 * ------------------------------------------------------------------------- */

void des_automaton_init(des_automaton_t *G, const char *name);

uint16_t des_automaton_add_state(des_automaton_t *G);

void des_automaton_mark_state(des_automaton_t *G, uint16_t q);

int des_automaton_is_marked(const des_automaton_t *G, uint16_t q);

uint16_t des_automaton_add_event(des_automaton_t *G, const char *label,
                                  uint8_t flags);

int des_automaton_add_transition(des_automaton_t *G, uint16_t src,
                                  uint16_t ev, uint16_t dst);

uint16_t des_automaton_delta(const des_automaton_t *G, uint16_t q,
                              uint16_t ev);

uint16_t des_automaton_delta_star(const des_automaton_t *G, uint16_t q,
                                   const uint16_t *s, size_t len);

int des_automaton_is_legal_string(const des_automaton_t *G, uint16_t q,
                                   const uint16_t *s, size_t len);

uint16_t des_automaton_epsilon_closure(const des_automaton_t *G, uint16_t q,
                                        uint16_t *reachable, uint16_t capacity);

uint16_t des_automaton_active_events(const des_automaton_t *G, uint16_t q,
                                      uint16_t *active, uint16_t cap);

int des_automaton_is_active(const des_automaton_t *G, uint16_t q,
                             uint16_t ev);

uint16_t des_automaton_controllable_at(const des_automaton_t *G, uint16_t q,
                                        uint16_t *ctrl, uint16_t cap);

uint16_t des_automaton_uncontrollable_at(const des_automaton_t *G, uint16_t q,
                                          uint16_t *unctrl, uint16_t cap);

uint16_t des_automaton_reachable_states(const des_automaton_t *G,
                                         uint16_t *reachable, uint16_t cap);

uint16_t des_automaton_coreachable_states(const des_automaton_t *G,
                                           uint16_t *coreach, uint16_t cap);

void des_automaton_trim(const des_automaton_t *G, des_automaton_t *trim);

int des_automaton_is_nonblocking(const des_automaton_t *G);

int des_automaton_is_prefix_closed(const des_automaton_t *G);

uint16_t des_automaton_generate_language(const des_automaton_t *G,
    uint16_t *strings, uint8_t *lens, uint16_t max_strs, uint8_t maxlen);

uint16_t des_automaton_marked_language(const des_automaton_t *G,
    uint16_t *strings, uint8_t *lens, uint16_t max_strs, uint8_t maxlen);

void des_automaton_product(const des_automaton_t *G1,
                            const des_automaton_t *G2,
                            des_automaton_t *prod);

size_t des_automaton_natural_projection(const uint16_t *s, size_t len,
    const des_automaton_t *G, uint16_t *proj, size_t proj_cap);

void des_automaton_print_dot(const des_automaton_t *G, FILE *fp);

void des_automaton_print_table(const des_automaton_t *G, FILE *fp);

uint16_t des_automaton_count_deadlocks(const des_automaton_t *G);

uint16_t des_automaton_count_livelocks(const des_automaton_t *G);

void des_automaton_destroy(des_automaton_t *G);

#endif /* DES_AUTOMATON_H */
