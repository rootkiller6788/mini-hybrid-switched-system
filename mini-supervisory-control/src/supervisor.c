/**
 * @file supervisor.c
 * @brief Supervisor Implementation - Ramadge-Wonham Framework
 *
 * Implements supervisor decision logic: S: L(G) -> Gamma
 * where Gamma = { gamma subset Sigma | Sigma_u subset gamma }
 *
 * Reference: Ramadge & Wonham, SIAM J. Control Optim. 25(1), 1987
 * Reference: Wonham & Cai (2019)
 */

#include "supervisor.h"
#include "controllability.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* ===================================================================
 * L1: Supervisor initialization
 * =================================================================== */

void supervisor_init(supervisor_t *S, des_automaton_t *plant,
                      des_automaton_t *spec, sup_type_t type,
                      const char *name)
{
    if (!S) return;
    memset(S, 0, sizeof(*S));
    S->plant = plant;
    S->spec  = spec;
    S->type  = type;
    S->current_state = 0;
    if (name) {
        strncpy(S->name, name, sizeof(S->name) - 1);
        S->name[sizeof(S->name) - 1] = '\0';
    } else {
        strcpy(S->name, "unnamed_sup");
    }
    S->npatterns = 0;
    S->hist_len = 0;
    S->nmodules = 0;
    S->events_observed = 0;
    S->events_disabled = 0;
    S->events_enabled = 0;
    S->violations = 0;
    if (plant) S->current_state = plant->q0;
}

/* ===================================================================
 * L2: Core supervisor operations
 * =================================================================== */

void supervisor_set_pattern(supervisor_t *S, uint16_t state,
                             uint32_t enabled)
{
    if (!S || S->npatterns >= SUP_MAX_PATTERNS) return;
    /* Check if pattern for this state already exists */
    for (uint16_t i = 0; i < S->npatterns; i++) {
        if (S->patterns[i].state == state) {
            S->patterns[i].enabled_mask = enabled;
            return;
        }
    }
    S->patterns[S->npatterns].state = state;
    S->patterns[S->npatterns].enabled_mask = enabled;
    S->npatterns++;
}

uint32_t supervisor_current_pattern(const supervisor_t *S, uint16_t state)
{
    if (!S) return 0;
    for (uint16_t i = 0; i < S->npatterns; i++) {
        if (S->patterns[i].state == state) {
            return S->patterns[i].enabled_mask;
        }
    }
    /* Default: enable all controllable events */
    return 0xFFFFFFFF;
}

int supervisor_is_enabled(const supervisor_t *S, uint16_t ev)
{
    if (!S || !S->plant) return 1; /* No plant = everything enabled */

    /* Uncontrollable events are always enabled (by definition) */
    if (ev < S->plant->nevents &&
        (S->plant->events[ev].flags & DES_EVENT_UNCONTROLLABLE)) {
        return 1;
    }

    /* Check control pattern for current state */
    uint32_t pattern = supervisor_current_pattern(S, S->current_state);
    return (pattern & (1u << ev)) ? 1 : 0;
}

int supervisor_observe(supervisor_t *S, uint16_t ev)
{
    if (!S) return 1;
    S->events_observed++;

    /* Record in history */
    if (S->hist_len < SUP_MAX_HISTORY) {
        S->history[S->hist_len++] = ev;
    }

    /* Check if this event is enabled by the supervisor.
     * MUST check BEFORE updating state because the control pattern
     * applies to the state where the event occurs. */
    if (!supervisor_is_enabled(S, ev)) {
        S->events_disabled++;
        S->violations++;
        return 0;
    }

    /* Check if event is legal in plant */
    if (S->plant) {
        if (!des_automaton_is_active(S->plant, S->current_state, ev)) {
            S->violations++;
            return 0;
        }
        /* Update current state */
        uint16_t nxt = des_automaton_delta(S->plant, S->current_state, ev);
        if (nxt == DES_UNDEF_STATE) {
            S->violations++;
            return 0;
        }
        S->current_state = nxt;
    }

    /* Check specification (if any) */
    if (S->spec) {
        if (!des_automaton_is_active(S->spec, S->current_state, ev)) {
            /* Spec tracking maintained for completeness */
        }
    }

    S->events_enabled++;
    return 1;
}

uint16_t supervisor_enabled_set(const supervisor_t *S, uint16_t *enabled,
                                 uint16_t cap)
{
    if (!S || !S->plant) return 0;
    uint16_t count = 0;
    uint16_t active[DES_MAX_EVENTS];
    uint16_t n_active = des_automaton_active_events(S->plant,
        S->current_state, active, DES_MAX_EVENTS);
    for (uint16_t i = 0; i < n_active && count < cap; i++) {
        if (supervisor_is_enabled(S, active[i])) {
            if (enabled) enabled[count] = active[i];
            count++;
        }
    }
    return count;
}

uint16_t supervisor_disabled_set(const supervisor_t *S, uint16_t *disabled,
                                  uint16_t cap)
{
    if (!S || !S->plant) return 0;
    uint16_t count = 0;
    uint16_t active[DES_MAX_EVENTS];
    uint16_t n_active = des_automaton_active_events(S->plant,
        S->current_state, active, DES_MAX_EVENTS);
    for (uint16_t i = 0; i < n_active && count < cap; i++) {
        /* Only controllable events can be disabled */
        if (!(S->plant->events[active[i]].flags & DES_EVENT_CONTROLLABLE))
            continue;
        if (!supervisor_is_enabled(S, active[i])) {
            if (disabled) disabled[count] = active[i];
            count++;
        }
    }
    return count;
}

int supervisor_is_proper(const supervisor_t *S)
{
    if (!S || !S->plant) return 0;
    /* Check every pattern: no uncontrollable event should be disabled */
    for (uint16_t i = 0; i < S->npatterns; i++) {
        uint32_t mask = S->patterns[i].enabled_mask;
        for (uint16_t ev = 0; ev < S->plant->nevents; ev++) {
            if ((S->plant->events[ev].flags & DES_EVENT_UNCONTROLLABLE) &&
                !(mask & (1u << ev))) {
                return 0; /* Uncontrollable event disabled! */
            }
        }
    }
    return 1;
}

void supervisor_closed_loop(const supervisor_t *S, des_automaton_t *closed_loop)
{
    if (!S || !S->plant || !closed_loop) return;

    char cl_name[128];
    snprintf(cl_name, sizeof(cl_name), "%s/%s", S->name, S->plant->name);
    des_automaton_init(closed_loop, cl_name);

    /* Copy all events from plant */
    for (uint16_t ev = 0; ev < S->plant->nevents; ev++) {
        des_automaton_add_event(closed_loop, S->plant->events[ev].label,
                                S->plant->events[ev].flags);
    }

    /* Copy all states with factory mappings */
    int16_t old2new[DES_MAX_STATES];
    for (uint16_t q = 0; q < S->plant->nstates; q++) {
        old2new[q] = (int16_t)des_automaton_add_state(closed_loop);
    }
    closed_loop->q0 = (uint16_t)old2new[S->plant->q0];

    /* Copy marked states */
    for (uint16_t i = 0; i < S->plant->n_marked; i++) {
        des_automaton_mark_state(closed_loop, (uint16_t)old2new[S->plant->marked[i]]);
    }

    /* Copy only transitions that the supervisor enables */
    for (uint16_t i = 0; i < S->plant->ntrans; i++) {
        uint16_t src = S->plant->transitions[i].src;
        uint16_t ev  = S->plant->transitions[i].ev;
        uint16_t dst = S->plant->transitions[i].dst;

        /* Check if supervisor enables this event at state src */
        uint32_t pattern = supervisor_current_pattern(S, src);
        if (pattern & (1u << ev)) {
            des_automaton_add_transition(closed_loop,
                (uint16_t)old2new[src], ev, (uint16_t)old2new[dst]);
        }
    }
}

/* ===================================================================
 * Modular supervision
 * =================================================================== */

int supervisor_add_module(supervisor_t *S, supervisor_t *module_sup)
{
    if (!S || !module_sup || S->nmodules >= SUP_MAX_MODULES) return 0;
    S->modules[S->nmodules++] = module_sup;
    return 1;
}

uint32_t supervisor_modular_pattern(const supervisor_t *S, uint16_t state)
{
    if (!S) return 0;
    if (S->nmodules == 0) return 0xFFFFFFFF;

    /* Intersection of all module decisions */
    uint32_t combined = 0xFFFFFFFF;
    for (uint16_t i = 0; i < S->nmodules; i++) {
        uint32_t mod_pattern = supervisor_current_pattern(S->modules[i], state);
        combined &= mod_pattern;
    }
    return combined;
}

void supervisor_reset(supervisor_t *S)
{
    if (!S) return;
    S->hist_len = 0;
    S->current_state = S->plant ? S->plant->q0 : 0;
    S->events_observed = 0;
    S->events_disabled = 0;
    S->events_enabled = 0;
    S->violations = 0;
}

void supervisor_print_status(const supervisor_t *S, FILE *fp)
{
    if (!S) return;
    if (!fp) fp = stdout;
    fprintf(fp, "=== Supervisor: %s ===\n", S->name);
    fprintf(fp, "Type: %d, Plant: %s\n", S->type,
            S->plant ? S->plant->name : "(none)");
    fprintf(fp, "Spec: %s\n", S->spec ? S->spec->name : "(none)");
    fprintf(fp, "Patterns: %u, Modules: %u\n", S->npatterns, S->nmodules);
    fprintf(fp, "Current state: %u, History len: %u\n",
            S->current_state, S->hist_len);
    fprintf(fp, "Stats: observed=%u, enabled=%u, disabled=%u, violations=%u\n",
            S->events_observed, S->events_enabled,
            S->events_disabled, S->violations);
    fprintf(fp, "Proper: %s\n", supervisor_is_proper(S) ? "yes" : "NO!");
    fprintf(fp, "Control patterns:\n");
    for (uint16_t i = 0; i < S->npatterns; i++) {
        fprintf(fp, "  state %u: mask=0x%08x\n",
                S->patterns[i].state, S->patterns[i].enabled_mask);
    }
}

void supervisor_destroy(supervisor_t *S)
{
    if (!S) return;
    memset(S, 0, sizeof(*S));
}