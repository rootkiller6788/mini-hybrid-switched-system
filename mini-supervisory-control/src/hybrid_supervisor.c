/**
 * @file hybrid_supervisor.c
 * @brief Hybrid Switched System Supervisor Implementation
 *
 * Bridges DES supervisory control (Ramadge-Wonham) with hybrid automata
 * for switched systems. The supervisor:
 *   1. Monitors continuous state evolution within each mode
 *   2. Guards mode transitions per safety specifications
 *   3. Uses DES abstraction for high-level decision making
 *   4. Enforces safety constraints: no unsafe mode reachable
 *
 * Applications: power electronics, flight control mode switching,
 *   manufacturing reconfiguration, autonomous vehicle mode management,
 *   building HVAC zone control (smart grid switching).
 *
 * Reference: Koutsoukos et al. (2000) "Supervisory control of hybrid systems"
 * Reference: Tabuada (2009) "Verification and control of hybrid systems"
 * Reference: Lygeros, Johansson et al. (2012)
 */

#include "hybrid_supervisor.h"
#include "controllability.h"
#include "synthesis.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* ===================================================================
 * L1: Hybrid automaton initialization
 * =================================================================== */

void hybrid_automaton_init(hybrid_automaton_t *H, const char *name)
{
    if (!H) return;
    memset(H, 0, sizeof(*H));

    if (name) {
        strncpy(H->name, name, sizeof(H->name) - 1);
        H->name[sizeof(H->name) - 1] = '\0';
    } else {
        strcpy(H->name, "unnamed_hybrid");
    }

    des_automaton_init(&H->des, name);
    H->nmodes = 0;
    H->ntrans = 0;
    H->current_mode = 0;
    H->t = 0.0;
    memset(&H->x, 0, sizeof(H->x));
    H->x.dim = 0;
}

uint16_t hybrid_add_mode(hybrid_automaton_t *H, const char *mode_name,
                          const hyb_ode_t *dynamics)
{
    if (!H || H->nmodes >= HYB_MAX_MODES) return HYB_MAX_MODES;

    uint16_t idx = H->nmodes;
    H->modes[idx].id = idx;
    H->modes[idx].name = mode_name;
    H->modes[idx].is_safe = 1;
    H->modes[idx].is_target = 0;
    H->modes[idx].n_invariants = 0;

    if (dynamics) {
        memcpy(&H->modes[idx].dynamics, dynamics, sizeof(hyb_ode_t));
    } else {
        memset(&H->modes[idx].dynamics, 0, sizeof(hyb_ode_t));
    }

    /* Also add a state in the DES automaton for this mode */
    des_automaton_add_state(&H->des);

    H->nmodes++;
    return idx;
}

void hybrid_set_invariant(hybrid_automaton_t *H, uint16_t mode,
                           const hyb_invariant_t *inv)
{
    if (!H || mode >= H->nmodes) return;
    if (!inv) return;

    if (H->modes[mode].n_invariants < HYB_MAX_INVARIANTS) {
        /* Store invariant: we store the pointer to a static or
         * user-managed invariant structure */
        memcpy(&H->modes[mode].invariant, inv, sizeof(hyb_invariant_t));
        H->modes[mode].n_invariants = 1;
    }
}

int hybrid_add_transition(hybrid_automaton_t *H, uint16_t src, uint16_t dst,
                           uint16_t event, hyb_switch_type_t sw_type,
                           const hyb_guard_t *guard,
                           const double reset_R[HYB_DIM][HYB_DIM],
                           const double *reset_r, uint8_t reset_dim)
{
    if (!H || src >= H->nmodes || dst >= H->nmodes ||
        H->ntrans >= HYB_MAX_GUARDS) return 0;

    uint16_t idx = H->ntrans;
    H->transitions[idx].src_mode = src;
    H->transitions[idx].dst_mode = dst;
    H->transitions[idx].event = event;
    H->transitions[idx].sw_type = sw_type;
    H->transitions[idx].reset_dim = reset_dim;

    if (guard) {
        memcpy(&H->transitions[idx].guard, guard, sizeof(hyb_guard_t));
    } else {
        memset(&H->transitions[idx].guard, 0, sizeof(hyb_guard_t));
    }

    if (reset_R) {
        memcpy(H->transitions[idx].reset_R, reset_R,
               HYB_DIM * HYB_DIM * sizeof(double));
    } else {
        memset(H->transitions[idx].reset_R, 0,
               HYB_DIM * HYB_DIM * sizeof(double));
        /* Default: identity reset */
        for (uint8_t d = 0; d < reset_dim && d < HYB_DIM; d++) {
            H->transitions[idx].reset_R[d][d] = 1.0;
        }
    }

    if (reset_r) {
        memcpy(H->transitions[idx].reset_r, reset_r,
               HYB_DIM * sizeof(double));
    } else {
        memset(H->transitions[idx].reset_r, 0, HYB_DIM * sizeof(double));
    }

    /* Also add transition in the DES layer */
    des_automaton_add_transition(&H->des, src, event, dst);

    H->ntrans++;
    return 1;
}

/* ===================================================================
 * L2: Hybrid operations - guards, invariants, continuous step
 * =================================================================== */

int hybrid_guard_satisfied(const hybrid_automaton_t *H, uint16_t trans_idx)
{
    if (!H || trans_idx >= H->ntrans) return 0;

    const hyb_guard_t *g = &H->transitions[trans_idx].guard;
    if (g->dim == 0) return 1; /* No guard = always satisfied */

    /* Linear inequality: A^T x + b <= 0 */
    double val = g->b;
    for (uint8_t d = 0; d < g->dim && d < H->x.dim; d++) {
        val += g->A[d] * H->x.data[d];
    }
    return (val <= 1e-10) ? 1 : 0;
}

int hybrid_invariant_satisfied(const hybrid_automaton_t *H, uint16_t mode)
{
    if (!H || mode >= H->nmodes) return 0;

    const hyb_invariant_t *inv = &H->modes[mode].invariant;
    if (inv->dim == 0) return 1; /* No invariant = always satisfied */

    /* Linear inequality: A^T x + b <= 0 */
    double val = inv->b;
    for (uint8_t d = 0; d < inv->dim && d < H->x.dim; d++) {
        val += inv->A[d] * H->x.data[d];
    }
    return (val <= 1e-10) ? 1 : 0;
}

void hybrid_step_continuous(hybrid_automaton_t *H, double dt,
                             const double *u, uint8_t u_dim)
{
    if (!H || H->current_mode >= H->nmodes) return;
    if (dt <= 0.0) return;

    const hyb_ode_t *ode = &H->modes[H->current_mode].dynamics;
    uint8_t dim = ode->dim;
    if (dim == 0 || dim > HYB_DIM) return;

    /* Forward Euler integration: x_{k+1} = x_k + dt * (A*x_k + B*u) */
    double dx[HYB_DIM];
    memset(dx, 0, sizeof(dx));

    /* dx = A*x */
    for (uint8_t i = 0; i < dim; i++) {
        for (uint8_t j = 0; j < dim; j++) {
            dx[i] += ode->A[i][j] * H->x.data[j];
        }
    }

    /* dx += B*u + c */
    if (u && u_dim > 0) {
        for (uint8_t i = 0; i < dim; i++) {
            for (uint8_t j = 0; j < u_dim && j < HYB_DIM; j++) {
                dx[i] += ode->B[i][j] * u[j];
            }
        }
    }
    for (uint8_t i = 0; i < dim; i++) {
        dx[i] += ode->c[i];
    }

    /* Update state */
    for (uint8_t i = 0; i < dim; i++) {
        H->x.data[i] += dt * dx[i];
    }

    H->x.dim = dim;
    H->t += dt;
}

int hybrid_execute_event(hybrid_automaton_t *H, uint16_t ev)
{
    if (!H || H->current_mode >= H->nmodes) return 0;

    /* Find a transition from current mode on this event */
    for (uint16_t i = 0; i < H->ntrans; i++) {
        const hyb_transition_t *tr = &H->transitions[i];
        if (tr->src_mode != H->current_mode) continue;
        if (tr->event != ev) continue;

        /* Check guard */
        if (!hybrid_guard_satisfied(H, i)) continue;

        /* Execute transition */
        uint16_t prev_mode = H->current_mode;
        H->current_mode = tr->dst_mode;

        /* Apply reset map */
        hybrid_reset_state(H, i);

        /* Also update DES state */
        uint16_t nxt_des = des_automaton_delta(&H->des, prev_mode, ev);
        if (nxt_des != DES_UNDEF_STATE) {
            /* DES tracks mode transitions */
        }

        return 1;
    }
    return 0; /* No enabled transition found */
}

void hybrid_reset_state(hybrid_automaton_t *H, uint16_t trans_idx)
{
    if (!H || trans_idx >= H->ntrans) return;

    const hyb_transition_t *tr = &H->transitions[trans_idx];
    uint8_t dim = tr->reset_dim;
    if (dim == 0 || dim > HYB_DIM) return;

    /* x' = R*x + r */
    double new_x[HYB_DIM];
    memset(new_x, 0, sizeof(new_x));

    for (uint8_t i = 0; i < dim; i++) {
        new_x[i] = tr->reset_r[i];
        for (uint8_t j = 0; j < dim; j++) {
            new_x[i] += tr->reset_R[i][j] * H->x.data[j];
        }
    }

    memcpy(H->x.data, new_x, HYB_DIM * sizeof(double));
    H->x.dim = dim;
}

uint16_t hybrid_enabled_transitions(const hybrid_automaton_t *H,
                                     uint16_t *enabled, uint16_t cap)
{
    if (!H || H->current_mode >= H->nmodes) return 0;

    uint16_t count = 0;
    for (uint16_t i = 0; i < H->ntrans && count < cap; i++) {
        if (H->transitions[i].src_mode != H->current_mode) continue;
        if (hybrid_guard_satisfied(H, i)) {
            if (enabled) enabled[count] = i;
            count++;
        }
    }
    return count;
}

/* ===================================================================
 * L5: Hybrid supervisor algorithms
 * =================================================================== */

int hybrid_supervisor_safe(const hybrid_automaton_t *H,
                            const supervisor_t *S)
{
    if (!H || !S) return 0;

    /* Check: current mode is safe AND all reachable modes via
     * uncontrollable events are safe AND supervisor enables at
     * least one safe transition from every reachable state */

    uint16_t mode = H->current_mode;
    if (!H->modes[mode].is_safe) return 0;

    /* Check uncontrollable transitions: must all lead to safe modes */
    for (uint16_t i = 0; i < H->ntrans; i++) {
        const hyb_transition_t *tr = &H->transitions[i];
        if (tr->src_mode != mode) continue;
        if (tr->sw_type == HYB_SWITCH_AUTONOMOUS) {
            /* Autonomous transitions cannot be prevented */
            if (hybrid_guard_satisfied(H, i) &&
                !H->modes[tr->dst_mode].is_safe) {
                return 0;
            }
        }
    }

    /* Check that supervisor enables at least one controlled transition
     * from each non-target safe mode (nonblocking condition) */
    if (!H->modes[mode].is_target) {
        /* Must have at least one enabled controlled transition */
        int has_enabled = 0;
        for (uint16_t i = 0; i < H->ntrans; i++) {
            const hyb_transition_t *tr = &H->transitions[i];
            if (tr->src_mode != mode) continue;
            if (tr->sw_type != HYB_SWITCH_CONTROLLED) continue;
            if (hybrid_guard_satisfied(H, i) &&
                supervisor_is_enabled(S, tr->event)) {
                has_enabled = 1;
                break;
            }
        }
        if (!has_enabled) {
            /* Check if invariant holds (can stay in this mode) */
            if (!hybrid_invariant_satisfied(H, mode)) return 0;
        }
    }

    return 1;
}

void hybrid_supervisor_synthesize(const hybrid_automaton_t *H,
                                   const des_automaton_t *safety_spec,
                                   supervisor_t *S)
{
    if (!H || !safety_spec || !S) return;

    /* Step 1: Abstract the hybrid automaton to a DES */
    des_automaton_t abs_G;
    hybrid_abstraction(H, &abs_G);

    /* Step 2: Check controllability of safety_spec wrt abs_G */
    controllability_result_t ctrl_result;
    (void)controllability_check(&abs_G, safety_spec, &ctrl_result);

    /* Step 3: Synthesize supervisor from the abstraction */
    des_automaton_t supC;
    controllability_compute_supremal(&abs_G, safety_spec, &supC);

    /* Step 4: Build the supervisor */
    supervisor_init(S, &abs_G, (des_automaton_t *)safety_spec,
                    SUP_TYPE_MONOLITHIC, H->name);

    /* For each state in supC, define the control pattern:
     * enable controllable events that stay within supC */
    for (uint16_t q = 0; q < supC.nstates; q++) {
        uint32_t enabled = 0;
        for (uint16_t ev = 0; ev < supC.nevents; ev++) {
            if (supC.delta[q][ev] != DES_UNDEF_STATE) {
                enabled |= (1u << ev);
            }
        }
        supervisor_set_pattern(S, q, enabled);
    }
}

void hybrid_abstraction(const hybrid_automaton_t *H,
                         des_automaton_t *abs_G)
{
    if (!H || !abs_G) return;

    /* Abstract the hybrid automaton to a pure DES:
     * - Modes become DES states
     * - Transitions become DES events
     * - Safety: modes marked as safe/unsafe
     * - Controllability: controlled switches = controllable events,
     *   autonomous switches = uncontrollable events */

    char abs_name[128];
    snprintf(abs_name, sizeof(abs_name), "abs(%s)", H->name);
    des_automaton_init(abs_G, abs_name);

    /* Create states for each mode */
    int16_t mode2state[HYB_MAX_MODES];
    for (uint16_t m = 0; m < H->nmodes; m++) {
        mode2state[m] = (int16_t)des_automaton_add_state(abs_G);
    }
    abs_G->q0 = (uint16_t)mode2state[0];

    /* Mark safe modes */
    for (uint16_t m = 0; m < H->nmodes; m++) {
        if (H->modes[m].is_safe) {
            des_automaton_mark_state(abs_G, (uint16_t)mode2state[m]);
        }
    }

    /* Create events for transitions (deduplicate by label pattern) */
    for (uint16_t i = 0; i < H->ntrans; i++) {
        const hyb_transition_t *tr = &H->transitions[i];
        /* Create event if not already existing */
        char ev_label[32];
        snprintf(ev_label, sizeof(ev_label), "sw_%u_%u",
                 tr->src_mode, tr->dst_mode);

        uint8_t flags = DES_EVENT_OBSERVABLE;
        if (tr->sw_type == HYB_SWITCH_CONTROLLED) {
            flags |= DES_EVENT_CONTROLLABLE;
        } else {
            flags |= DES_EVENT_UNCONTROLLABLE;
        }

        /* Check if event with this label already exists */
        int16_t ev_idx = -1;
        for (uint16_t e = 0; e < abs_G->nevents; e++) {
            if (abs_G->events[e].label &&
                strcmp(abs_G->events[e].label, ev_label) == 0) {
                ev_idx = (int16_t)e;
                break;
            }
        }
        if (ev_idx < 0) {
            ev_idx = (int16_t)des_automaton_add_event(abs_G, ev_label, flags);
        }

        /* Add transition */
        des_automaton_add_transition(abs_G,
            (uint16_t)mode2state[tr->src_mode],
            (uint16_t)ev_idx,
            (uint16_t)mode2state[tr->dst_mode]);
    }
}

/* ===================================================================
 * L5: Hybrid simulation
 * =================================================================== */

int hybrid_simulate(hybrid_automaton_t *H, supervisor_t *S,
                     double T_end, double dt,
                     hyb_state_t *trajectory, uint16_t max_steps,
                     uint16_t *n_steps)
{
    if (!H || !H->nmodes) return 0;

    uint16_t steps = 0;
    double t = H->t;

    /* Record initial state */
    if (trajectory && steps < max_steps) {
        memcpy(&trajectory[steps], &H->x, sizeof(hyb_state_t));
        steps++;
    }

    while (t < T_end && steps < max_steps) {
        uint16_t mode = H->current_mode;

        /* Check invariant */
        if (!hybrid_invariant_satisfied(H, mode)) {
            /* Must take an autonomous transition */
            int took_transition = 0;
            for (uint16_t i = 0; i < H->ntrans; i++) {
                if (H->transitions[i].src_mode != mode) continue;
                if (H->transitions[i].sw_type != HYB_SWITCH_AUTONOMOUS) continue;
                if (hybrid_guard_satisfied(H, i)) {
                    H->current_mode = H->transitions[i].dst_mode;
                    hybrid_reset_state(H, i);
                    took_transition = 1;
                    break;
                }
            }
            if (!took_transition) {
                /* Invariant violated, no autonomous transition available:
                 * system enters unsafe state */
                H->modes[mode].is_safe = 0;
                break;
            }
        }

        /* Check for supervisor-enabled controlled transitions */
        uint16_t enabled[HYB_MAX_GUARDS];
        uint16_t n_enabled = hybrid_enabled_transitions(H, enabled,
                                                         HYB_MAX_GUARDS);

        int took_controlled = 0;
        for (uint16_t i = 0; i < n_enabled && !took_controlled; i++) {
            uint16_t tr_idx = enabled[i];
            const hyb_transition_t *tr = &H->transitions[tr_idx];

            if (tr->sw_type != HYB_SWITCH_CONTROLLED) continue;
            if (!supervisor_is_enabled(S, tr->event)) continue;

            /* Take the first enabled controlled transition */
            supervisor_observe(S, tr->event);
            H->current_mode = tr->dst_mode;
            hybrid_reset_state(H, tr_idx);
            took_controlled = 1;
        }

        /* Continuous evolution */
        if (!took_controlled) {
            hybrid_step_continuous(H, dt, NULL, 0);
        }

        /* Record state */
        if (trajectory && steps < max_steps) {
            memcpy(&trajectory[steps], &H->x, sizeof(hyb_state_t));
            steps++;
        }

        t += dt;
    }

    if (n_steps) *n_steps = steps;
    return 1;
}

/* ===================================================================
 * L5: Reachability check
 * =================================================================== */

int hybrid_reachability_check(const hybrid_automaton_t *H,
                               const hyb_state_t *init_set,
                               const hyb_state_t *unsafe_set)
{
    if (!H || !init_set || !unsafe_set) return 0;

    /* Simple forward reachability:
     * Check if there exists a sequence of modes and transitions
     * such that the unsafe set is reachable from the initial set.
     *
     * This is a very simplified check; full reachability requires
     * tools like PHAVer, SpaceEx, or Flow*.
     *
     * Here we check: is there a mode that is both reachable via
     * DES transitions and marked as unsafe? */

    /* First check: direct DES reachability of unsafe modes */
    for (uint16_t m = 0; m < H->nmodes; m++) {
        if (!H->modes[m].is_safe) {
            /* Check if this mode is reachable from q0 in the DES */
            uint16_t reachable[DES_MAX_STATES];
            uint16_t n_reach = des_automaton_reachable_states(&H->des,
                reachable, DES_MAX_STATES);
            for (uint16_t i = 0; i < n_reach; i++) {
                if (reachable[i] == m) {
                    return 0; /* Unsafe mode is reachable! */
                }
            }
        }
    }

    /* Second check: can the continuous dynamics drive the state
     * into the unsafe region without changing mode?
     * (simplified: check if any unsafe mode is adjacent via
     *  uncontrollable transitions) */
    for (uint16_t m = 0; m < H->nmodes; m++) {
        if (H->modes[m].is_safe) {
            for (uint16_t i = 0; i < H->ntrans; i++) {
                if (H->transitions[i].src_mode == m &&
                    H->transitions[i].sw_type == HYB_SWITCH_AUTONOMOUS) {
                    uint16_t dst = H->transitions[i].dst_mode;
                    if (!H->modes[dst].is_safe) {
                        return 0; /* Unsafe mode reachable via autonomous transition */
                    }
                }
            }
        }
    }

    return 1; /* Safe */
}

/* ===================================================================
 * L7: Utility functions
 * =================================================================== */

void hybrid_print_mode(const hybrid_automaton_t *H, uint16_t mode, FILE *fp)
{
    if (!H || mode >= H->nmodes) return;
    if (!fp) fp = stdout;

    const hyb_mode_t *m = &H->modes[mode];
    fprintf(fp, "Mode %u: %s\n", mode, m->name ? m->name : "?");
    fprintf(fp, "  Safe: %s, Target: %s\n",
            m->is_safe ? "yes" : "no",
            m->is_target ? "yes" : "no");
    fprintf(fp, "  Dynamics: dim=%u\n", m->dynamics.dim);
    if (m->dynamics.dim > 0) {
        fprintf(fp, "  A matrix:\n");
        for (uint8_t i = 0; i < m->dynamics.dim; i++) {
            fprintf(fp, "    ");
            for (uint8_t j = 0; j < m->dynamics.dim; j++) {
                fprintf(fp, "%8.3f ", m->dynamics.A[i][j]);
            }
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "  Invariants: %u\n", m->n_invariants);
}

void hybrid_print_trajectory(const hyb_state_t *traj, uint16_t n_steps,
                              FILE *fp)
{
    if (!traj || !fp) return;

    fprintf(fp, "Trajectory (%u steps):\n", n_steps);
    for (uint16_t i = 0; i < n_steps; i++) {
        fprintf(fp, "  [%4u] ", i);
        for (uint8_t d = 0; d < traj[i].dim && d < HYB_DIM; d++) {
            fprintf(fp, "%10.6f ", traj[i].data[d]);
        }
        fprintf(fp, "\n");
    }
}

void hybrid_automaton_destroy(hybrid_automaton_t *H)
{
    if (!H) return;
    memset(H, 0, sizeof(*H));
}