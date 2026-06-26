/**
 * @file des_automaton.c
 * @brief DES Automaton Implementation - Core Operations
 *
 * Implements the deterministic finite automaton model for Discrete Event Systems
 * in the Ramadge-Wonham supervisory control framework.
 *
 * Reference: Ramadge & Wonham (1989), Cassandras & Lafortune (2008)
 */

#include "des_automaton.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* ===================================================================
 * L1: Core Definitions - Init, State/Event/Transition Management
 * =================================================================== */

void des_automaton_init(des_automaton_t *G, const char *name)
{
    if (!G) return;
    memset(G, 0, sizeof(*G));
    if (name) {
        strncpy(G->name, name, sizeof(G->name) - 1);
        G->name[sizeof(G->name) - 1] = '\0';
    } else {
        strcpy(G->name, "unnamed");
    }
    G->q0 = 0;
    G->nstates = 1;
    G->nevents = 0;
    G->ntrans = 0;
    G->n_controllable = 0;
    G->n_observable = 0;
    G->n_marked = 0;
    for (uint16_t i = 0; i < DES_MAX_STATES; i++) {
        for (uint16_t j = 0; j < DES_MAX_EVENTS; j++) {
            G->delta[i][j] = DES_UNDEF_STATE;
        }
    }
}

uint16_t des_automaton_add_state(des_automaton_t *G)
{
    if (!G || G->nstates >= DES_MAX_STATES) return DES_UNDEF_STATE;
    uint16_t idx = G->nstates;
    G->nstates++;
    for (uint16_t j = 0; j < DES_MAX_EVENTS; j++) {
        G->delta[idx][j] = DES_UNDEF_STATE;
    }
    return idx;
}

void des_automaton_mark_state(des_automaton_t *G, uint16_t q)
{
    if (!G || q >= G->nstates) return;
    if (G->n_marked >= DES_MAX_STATES) return;
    for (uint16_t i = 0; i < G->n_marked; i++) {
        if (G->marked[i] == q) return;
    }
    G->marked[G->n_marked++] = q;
}

int des_automaton_is_marked(const des_automaton_t *G, uint16_t q)
{
    if (!G || q >= G->nstates) return 0;
    for (uint16_t i = 0; i < G->n_marked; i++) {
        if (G->marked[i] == q) return 1;
    }
    return 0;
}

uint16_t des_automaton_add_event(des_automaton_t *G, const char *label,
                                  uint8_t flags)
{
    if (!G || G->nevents >= DES_MAX_EVENTS) return DES_MAX_EVENTS;
    uint16_t idx = G->nevents;
    G->events[idx].idx   = idx;
    G->events[idx].label = label;
    G->events[idx].flags = flags;
    G->events[idx].cost  = 1.0;
    G->nevents++;
    if (flags & DES_EVENT_CONTROLLABLE) G->n_controllable++;
    if (flags & DES_EVENT_OBSERVABLE)   G->n_observable++;
    return idx;
}

int des_automaton_add_transition(des_automaton_t *G, uint16_t src,
                                  uint16_t ev, uint16_t dst)
{
    if (!G || src >= G->nstates || ev >= G->nevents
        || dst >= G->nstates) return 0;
    if (G->delta[src][ev] != DES_UNDEF_STATE) return 0;
    if (G->ntrans >= DES_MAX_TRANS) return 0;
    G->delta[src][ev] = dst;
    G->transitions[G->ntrans].src = src;
    G->transitions[G->ntrans].ev  = ev;
    G->transitions[G->ntrans].dst = dst;
    G->ntrans++;
    return 1;
}

uint16_t des_automaton_delta(const des_automaton_t *G, uint16_t q,
                              uint16_t ev)
{
    if (!G || q >= G->nstates || ev >= G->nevents)
        return DES_UNDEF_STATE;
    return G->delta[q][ev];
}

uint16_t des_automaton_delta_star(const des_automaton_t *G, uint16_t q,
                                   const uint16_t *s, size_t len)
{
    if (!G || q >= G->nstates || !s) return DES_UNDEF_STATE;
    uint16_t current = q;
    for (size_t i = 0; i < len; i++) {
        if (s[i] >= G->nevents) return DES_UNDEF_STATE;
        current = G->delta[current][s[i]];
        if (current == DES_UNDEF_STATE) return DES_UNDEF_STATE;
    }
    return current;
}

int des_automaton_is_legal_string(const des_automaton_t *G, uint16_t q,
                                   const uint16_t *s, size_t len)
{
    uint16_t end_state = des_automaton_delta_star(G, q, s, len);
    return (end_state != DES_UNDEF_STATE) ? 1 : 0;
}

/* ===================================================================
 * L2: Active events, reachability, epsilon-closure
 * =================================================================== */

uint16_t des_automaton_active_events(const des_automaton_t *G, uint16_t q,
                                      uint16_t *active, uint16_t cap)
{
    if (!G || q >= G->nstates) return 0;
    uint16_t count = 0;
    for (uint16_t ev = 0; ev < G->nevents && count < cap; ev++) {
        if (G->delta[q][ev] != DES_UNDEF_STATE) {
            if (active) active[count] = ev;
            count++;
        }
    }
    return count;
}

int des_automaton_is_active(const des_automaton_t *G, uint16_t q,
                             uint16_t ev)
{
    if (!G || q >= G->nstates || ev >= G->nevents) return 0;
    return (G->delta[q][ev] != DES_UNDEF_STATE) ? 1 : 0;
}

uint16_t des_automaton_controllable_at(const des_automaton_t *G, uint16_t q,
                                        uint16_t *ctrl, uint16_t cap)
{
    if (!G || q >= G->nstates) return 0;
    uint16_t count = 0;
    for (uint16_t ev = 0; ev < G->nevents && count < cap; ev++) {
        if (G->delta[q][ev] != DES_UNDEF_STATE &&
            (G->events[ev].flags & DES_EVENT_CONTROLLABLE)) {
            if (ctrl) ctrl[count] = ev;
            count++;
        }
    }
    return count;
}

uint16_t des_automaton_uncontrollable_at(const des_automaton_t *G, uint16_t q,
                                          uint16_t *unctrl, uint16_t cap)
{
    if (!G || q >= G->nstates) return 0;
    uint16_t count = 0;
    for (uint16_t ev = 0; ev < G->nevents && count < cap; ev++) {
        if (G->delta[q][ev] != DES_UNDEF_STATE &&
            (G->events[ev].flags & DES_EVENT_UNCONTROLLABLE)) {
            if (unctrl) unctrl[count] = ev;
            count++;
        }
    }
    return count;
}

uint16_t des_automaton_epsilon_closure(const des_automaton_t *G, uint16_t q,
                                        uint16_t *reachable, uint16_t capacity)
{
    if (!G || q >= G->nstates || !reachable || capacity == 0) return 0;
    uint8_t visited[DES_MAX_STATES];
    uint16_t queue[DES_MAX_STATES];
    uint16_t front = 0, back = 0, count = 0;
    memset(visited, 0, sizeof(visited));
    queue[back++] = q;
    visited[q] = 1;
    reachable[count++] = q;
    while (front < back && count < capacity) {
        uint16_t cur = queue[front++];
        for (uint16_t ev = 0; ev < G->nevents; ev++) {
            if (!(G->events[ev].flags & DES_EVENT_UNOBSERVABLE)) continue;
            uint16_t nxt = G->delta[cur][ev];
            if (nxt != DES_UNDEF_STATE && !visited[nxt]) {
                visited[nxt] = 1;
                queue[back++] = nxt;
                reachable[count++] = nxt;
                if (count >= capacity) return count;
            }
        }
    }
    return count;
}

/* ===================================================================
 * L4: Fundamental Laws - Reachability, Trim, Nonblocking
 * =================================================================== */

uint16_t des_automaton_reachable_states(const des_automaton_t *G,
                                         uint16_t *reachable, uint16_t cap)
{
    if (!G) return 0;
    uint8_t visited[DES_MAX_STATES];
    uint16_t queue[DES_MAX_STATES];
    uint16_t front = 0, back = 0, count = 0;
    memset(visited, 0, sizeof(visited));
    queue[back++] = G->q0;
    visited[G->q0] = 1;
    while (front < back) {
        uint16_t cur = queue[front++];
        if (reachable && count < cap) {
            reachable[count] = cur;
        }
        count++;
        for (uint16_t ev = 0; ev < G->nevents; ev++) {
            uint16_t nxt = G->delta[cur][ev];
            if (nxt != DES_UNDEF_STATE && !visited[nxt]) {
                visited[nxt] = 1;
                queue[back++] = nxt;
            }
        }
    }
    return count;
}

uint16_t des_automaton_coreachable_states(const des_automaton_t *G,
                                           uint16_t *coreach, uint16_t cap)
{
    if (!G) return 0;
    uint8_t visiting[DES_MAX_STATES];
    uint16_t queue[DES_MAX_STATES];
    uint16_t front = 0, back = 0, count = 0;
    memset(visiting, 0, sizeof(visiting));
    /* BFS backward from marked states */
    for (uint16_t i = 0; i < G->n_marked; i++) {
        uint16_t mq = G->marked[i];
        if (!visiting[mq]) {
            visiting[mq] = 1;
            queue[back++] = mq;
        }
    }
    while (front < back) {
        uint16_t cur = queue[front++];
        if (coreach && count < cap) {
            coreach[count] = cur;
        }
        count++;
        for (uint16_t src = 0; src < G->nstates; src++) {
            if (visiting[src]) continue;
            for (uint16_t ev = 0; ev < G->nevents; ev++) {
                if (G->delta[src][ev] == cur) {
                    visiting[src] = 1;
                    queue[back++] = src;
                    break;
                }
            }
        }
    }
    return count;
}

void des_automaton_trim(const des_automaton_t *G, des_automaton_t *trim)
{
    if (!G || !trim) return;
    uint16_t reachable[DES_MAX_STATES];
    uint8_t  is_reachable[DES_MAX_STATES];
    uint16_t coreachable[DES_MAX_STATES];
    uint8_t  is_coreachable[DES_MAX_STATES];
    uint8_t  keep[DES_MAX_STATES];
    int16_t  old2new[DES_MAX_STATES];
    memset(is_reachable, 0, sizeof(is_reachable));
    memset(is_coreachable, 0, sizeof(is_coreachable));
    memset(keep, 0, sizeof(keep));

    uint16_t n_reach = des_automaton_reachable_states(G, reachable, DES_MAX_STATES);
    for (uint16_t i = 0; i < n_reach; i++) is_reachable[reachable[i]] = 1;

    uint16_t n_coreach = des_automaton_coreachable_states(G, coreachable, DES_MAX_STATES);
    for (uint16_t i = 0; i < n_coreach; i++) is_coreachable[coreachable[i]] = 1;

    for (uint16_t q = 0; q < G->nstates; q++) {
        old2new[q] = -1;
        if (is_reachable[q] && is_coreachable[q]) keep[q] = 1;
    }

    des_automaton_init(trim, G->name);
    for (uint16_t ev = 0; ev < G->nevents; ev++) {
        des_automaton_add_event(trim, G->events[ev].label, G->events[ev].flags);
    }
    /* Reuse pre-allocated state 0 for G->q0 if it's kept */
    if (keep[G->q0]) {
        old2new[G->q0] = 0;
    }
    for (uint16_t q = 0; q < G->nstates; q++) {
        if (q == G->q0) continue; /* Already handled */
        if (keep[q]) old2new[q] = (int16_t)des_automaton_add_state(trim);
    }
    trim->q0 = (uint16_t)old2new[G->q0];
    for (uint16_t i = 0; i < G->n_marked; i++) {
        uint16_t mq = G->marked[i];
        if (keep[mq]) des_automaton_mark_state(trim, (uint16_t)old2new[mq]);
    }
    for (uint16_t i = 0; i < G->ntrans; i++) {
        uint16_t src = G->transitions[i].src;
        uint16_t ev  = G->transitions[i].ev;
        uint16_t dst = G->transitions[i].dst;
        if (keep[src] && keep[dst]) {
            des_automaton_add_transition(trim, (uint16_t)old2new[src], ev, (uint16_t)old2new[dst]);
        }
    }
}

int des_automaton_is_nonblocking(const des_automaton_t *G)
{
    if (!G) return 0;
    uint16_t n_reach = des_automaton_reachable_states(G, NULL, 0);
    uint16_t n_coreach = des_automaton_coreachable_states(G, NULL, 0);
    return (n_reach == n_coreach) ? 1 : 0;
}

int des_automaton_is_prefix_closed(const des_automaton_t *G)
{
    if (!G) return 0;
    uint16_t reachable[DES_MAX_STATES];
    uint16_t n_reach = des_automaton_reachable_states(G, reachable, DES_MAX_STATES);
    for (uint16_t i = 0; i < G->ntrans; i++) {
        uint16_t src = G->transitions[i].src;
        int reachable_found = 0;
        for (uint16_t j = 0; j < n_reach; j++) {
            if (reachable[j] == src) { reachable_found = 1; break; }
        }
        if (!reachable_found) return 0;
    }
    return 1;
}

/* ===================================================================
 * L4: Language generation and product
 * =================================================================== */

static void gen_recursive(const des_automaton_t *G, uint16_t q,
    uint16_t *current, uint8_t cur_len, uint8_t maxlen,
    uint16_t *strings, uint8_t *lens, uint16_t *count, uint16_t max_strs)
{
    if (*count >= max_strs || cur_len >= maxlen) return;
    for (uint16_t ev = 0; ev < G->nevents; ev++) {
        uint16_t nxt = G->delta[q][ev];
        if (nxt == DES_UNDEF_STATE) continue;
        current[cur_len] = ev;
        uint16_t idx = *count;
        if (idx < max_strs) {
            for (uint8_t k = 0; k <= cur_len; k++) {
                strings[idx * DES_MAX_EVENTS + k] = current[k];
            }
            lens[idx] = cur_len + 1;
        }
        (*count)++;
        gen_recursive(G, nxt, current, cur_len + 1, maxlen, strings, lens, count, max_strs);
    }
}

uint16_t des_automaton_generate_language(const des_automaton_t *G,
    uint16_t *strings, uint8_t *lens, uint16_t max_strs, uint8_t maxlen)
{
    if (!G || !max_strs) return 0;
    uint16_t current[DES_MAX_EVENTS];
    uint16_t count = 0;
    if (count < max_strs) { lens[count] = 0; count++; }
    gen_recursive(G, G->q0, current, 0, maxlen, strings, lens, &count, max_strs);
    return count;
}

uint16_t des_automaton_marked_language(const des_automaton_t *G,
    uint16_t *strings, uint8_t *lens, uint16_t max_strs, uint8_t maxlen)
{
    if (!G || !max_strs) return 0;
    uint16_t all_strs[4096];
    uint8_t  all_lens[4096];
    uint16_t n_all = des_automaton_generate_language(G, all_strs, all_lens, 4096, maxlen);
    uint16_t count = 0;
    for (uint16_t i = 0; i < n_all && count < max_strs; i++) {
        uint16_t end_state;
        if (all_lens[i] == 0) {
            end_state = G->q0;
        } else {
            end_state = des_automaton_delta_star(G, G->q0, &all_strs[i * DES_MAX_EVENTS], all_lens[i]);
        }
        if (end_state != DES_UNDEF_STATE && des_automaton_is_marked(G, end_state)) {
            for (uint8_t k = 0; k < all_lens[i]; k++) {
                strings[count * DES_MAX_EVENTS + k] = all_strs[i * DES_MAX_EVENTS + k];
            }
            lens[count] = all_lens[i];
            count++;
        }
    }
    return count;
}

void des_automaton_product(const des_automaton_t *G1,
                            const des_automaton_t *G2,
                            des_automaton_t *prod)
{
    if (!G1 || !G2 || !prod) return;

    char prod_name[128];
    snprintf(prod_name, sizeof(prod_name), "%s||%s", G1->name, G2->name);
    des_automaton_init(prod, prod_name);

    /* Event mapping */
    int16_t ev_map1[DES_MAX_EVENTS];
    int16_t ev_map2[DES_MAX_EVENTS];
    for (uint16_t i = 0; i < DES_MAX_EVENTS; i++) { ev_map1[i] = -1; ev_map2[i] = -1; }

    for (uint16_t i = 0; i < G1->nevents; i++) {
        ev_map1[i] = (int16_t)des_automaton_add_event(prod, G1->events[i].label, G1->events[i].flags);
    }
    for (uint16_t j = 0; j < G2->nevents; j++) {
        int found = 0;
        for (uint16_t i = 0; i < G1->nevents; i++) {
            if (G1->events[i].label && G2->events[j].label &&
                strcmp(G1->events[i].label, G2->events[j].label) == 0) {
                ev_map2[j] = ev_map1[i];
                found = 1;
                break;
            }
        }
        if (!found) {
            ev_map2[j] = (int16_t)des_automaton_add_event(prod, G2->events[j].label, G2->events[j].flags);
        }
    }

    /* State pair mapping */
    int16_t pair_map[DES_MAX_STATES][DES_MAX_STATES];
    for (uint16_t i = 0; i < G1->nstates; i++)
        for (uint16_t j = 0; j < G2->nstates; j++)
            pair_map[i][j] = -1;

    uint16_t queue1[DES_MAX_STATES * DES_MAX_STATES];
    uint16_t queue2[DES_MAX_STATES * DES_MAX_STATES];
    uint16_t front = 0, back = 0;

    queue1[back] = G1->q0;
    queue2[back] = G2->q0;
    back++;

    uint16_t new_q = des_automaton_add_state(prod);
    pair_map[G1->q0][G2->q0] = (int16_t)new_q;
    prod->q0 = new_q;

    if (des_automaton_is_marked(G1, G1->q0) && des_automaton_is_marked(G2, G2->q0)) {
        des_automaton_mark_state(prod, new_q);
    }

    while (front < back) {
        uint16_t q1 = queue1[front];
        uint16_t q2 = queue2[front];
        front++;
        int16_t prod_src = pair_map[q1][q2];
        if (prod_src < 0) continue;

        for (uint16_t ev = 0; ev < prod->nevents; ev++) {
            uint16_t dst1 = q1, dst2 = q2;
            int valid = 0;
            for (uint16_t i = 0; i < G1->nevents; i++) {
                if (ev_map1[i] == (int16_t)ev && G1->delta[q1][i] != DES_UNDEF_STATE) {
                    dst1 = G1->delta[q1][i]; valid |= 1;
                }
            }
            for (uint16_t j = 0; j < G2->nevents; j++) {
                if (ev_map2[j] == (int16_t)ev && G2->delta[q2][j] != DES_UNDEF_STATE) {
                    dst2 = G2->delta[q2][j]; valid |= 2;
                }
            }
            int is_shared = 0;
            for (uint16_t i = 0; i < G1->nevents && !is_shared; i++) {
                if (ev_map1[i] == (int16_t)ev) {
                    for (uint16_t j = 0; j < G2->nevents; j++) {
                        if (ev_map2[j] == (int16_t)ev) { is_shared = 1; break; }
                    }
                }
            }
            if (is_shared && valid != 3) continue;
            if (!is_shared && valid == 0) continue;

            if (pair_map[dst1][dst2] < 0) {
                uint16_t new_s = des_automaton_add_state(prod);
                pair_map[dst1][dst2] = (int16_t)new_s;
                if (des_automaton_is_marked(G1, dst1) && des_automaton_is_marked(G2, dst2)) {
                    des_automaton_mark_state(prod, new_s);
                }
                queue1[back] = dst1;
                queue2[back] = dst2;
                back++;
            }
            des_automaton_add_transition(prod, (uint16_t)prod_src, ev, (uint16_t)pair_map[dst1][dst2]);
        }
    }
}

size_t des_automaton_natural_projection(const uint16_t *s, size_t len,
    const des_automaton_t *G, uint16_t *proj, size_t proj_cap)
{
    if (!s || !G || !proj) return 0;
    size_t out_len = 0;
    for (size_t i = 0; i < len && out_len < proj_cap; i++) {
        uint16_t ev = s[i];
        if (ev < G->nevents && (G->events[ev].flags & DES_EVENT_OBSERVABLE)) {
            proj[out_len++] = ev;
        }
    }
    return out_len;
}

/* ===================================================================
 * Utility: print, deadlock, livelock
 * =================================================================== */

void des_automaton_print_dot(const des_automaton_t *G, FILE *fp)
{
    if (!G || !fp) return;
    fprintf(fp, "digraph %s {\n", G->name);
    fprintf(fp, "  rankdir=LR;\n  node [shape=circle];\n");
    for (uint16_t i = 0; i < G->n_marked; i++)
        fprintf(fp, "  %d [shape=doublecircle];\n", G->marked[i]);
    fprintf(fp, "  start [shape=point];\n  start -> %d;\n", G->q0);
    for (uint16_t i = 0; i < G->ntrans; i++) {
        const char *lbl = G->events[G->transitions[i].ev].label;
        const char *style = "";
        if (G->events[G->transitions[i].ev].flags & DES_EVENT_UNCONTROLLABLE)
            style = " [style=dashed]";
        fprintf(fp, "  %d -> %d [label=\"%s\"%s];\n",
                G->transitions[i].src, G->transitions[i].dst,
                lbl ? lbl : "?", style);
    }
    fprintf(fp, "}\n");
}

void des_automaton_print_table(const des_automaton_t *G, FILE *fp)
{
    if (!G || !fp) return;
    fprintf(fp, "Automaton: %s\n", G->name);
    fprintf(fp, "States: %u (q0=%u, marked=%u)\n", G->nstates, G->q0, G->n_marked);
    fprintf(fp, "Events: %u (ctrl=%u, obs=%u)\n", G->nevents, G->n_controllable, G->n_observable);
    fprintf(fp, "Transitions: %u\n\n", G->ntrans);
    fprintf(fp, "Marked: ");
    for (uint16_t i = 0; i < G->n_marked; i++) fprintf(fp, "%u ", G->marked[i]);
    fprintf(fp, "\n\nEvents:\n");
    for (uint16_t i = 0; i < G->nevents; i++)
        fprintf(fp, "  [%u] %s (flags=0x%02x)\n", i, G->events[i].label ? G->events[i].label : "?", G->events[i].flags);
    fprintf(fp, "\nDelta:\n     ");
    for (uint16_t ev = 0; ev < G->nevents; ev++) fprintf(fp, "%5u", ev);
    fprintf(fp, "\n");
    for (uint16_t q = 0; q < G->nstates; q++) {
        fprintf(fp, " %3u:", q);
        for (uint16_t ev = 0; ev < G->nevents; ev++) {
            if (G->delta[q][ev] == DES_UNDEF_STATE) fprintf(fp, "    -");
            else fprintf(fp, "%5u", G->delta[q][ev]);
        }
        fprintf(fp, "\n");
    }
}

uint16_t des_automaton_count_deadlocks(const des_automaton_t *G)
{
    if (!G) return 0;
    uint16_t deadlocks = 0;
    for (uint16_t q = 0; q < G->nstates; q++) {
        int has_out = 0;
        for (uint16_t ev = 0; ev < G->nevents; ev++) {
            if (G->delta[q][ev] != DES_UNDEF_STATE) { has_out = 1; break; }
        }
        if (!has_out && !des_automaton_is_marked(G, q)) deadlocks++;
    }
    return deadlocks;
}

uint16_t des_automaton_count_livelocks(const des_automaton_t *G)
{
    if (!G) return 0;
    /* Tarjan SCC-based livelock detection */
    uint16_t index[DES_MAX_STATES], lowlink[DES_MAX_STATES];
    int      onstack[DES_MAX_STATES];
    uint16_t stack[DES_MAX_STATES], stack_ptr = 0, cur_index = 0;
    uint16_t scc_id[DES_MAX_STATES], scc_count = 0;
    uint16_t scc_size[DES_MAX_STATES];
    int      scc_has_marked[DES_MAX_STATES];
    typedef struct { uint16_t q; uint16_t ev; uint8_t phase; } dfs_frame_t;
    dfs_frame_t dfs[DES_MAX_STATES * 2];
    uint16_t dfs_ptr = 0;

    memset(onstack, 0, sizeof(onstack));
    memset(scc_size, 0, sizeof(scc_size));
    memset(scc_has_marked, 0, sizeof(scc_has_marked));
    for (uint16_t i = 0; i < G->nstates; i++) {
        index[i] = DES_UNDEF_STATE;
        scc_id[i] = DES_UNDEF_STATE;
    }

    for (uint16_t start = 0; start < G->nstates; start++) {
        if (index[start] != DES_UNDEF_STATE) continue;
        dfs[0].q = start; dfs[0].ev = 0; dfs[0].phase = 0; dfs_ptr = 1;
        while (dfs_ptr > 0) {
            dfs_frame_t *cur = &dfs[dfs_ptr - 1];
            if (cur->phase == 0) {
                index[cur->q] = cur_index;
                lowlink[cur->q] = cur_index;
                cur_index++;
                stack[stack_ptr++] = cur->q;
                onstack[cur->q] = 1;
                cur->phase = 1;
            }
            int pushed = 0;
            while (cur->ev < G->nevents) {
                uint16_t nxt = G->delta[cur->q][cur->ev];
                cur->ev++;
                if (nxt == DES_UNDEF_STATE) continue;
                if (index[nxt] == DES_UNDEF_STATE) {
                    dfs[dfs_ptr].q = nxt; dfs[dfs_ptr].ev = 0;
                    dfs[dfs_ptr].phase = 0; dfs_ptr++; pushed = 1; break;
                } else if (onstack[nxt]) {
                    if (index[nxt] < lowlink[cur->q]) lowlink[cur->q] = index[nxt];
                }
            }
            if (!pushed) {
                dfs_ptr--;
                if (lowlink[cur->q] == index[cur->q]) {
                    uint16_t w; uint16_t sz = 0;
                    do {
                        w = stack[--stack_ptr];
                        onstack[w] = 0;
                        scc_id[w] = scc_count;
                        sz++;
                        if (des_automaton_is_marked(G, w)) scc_has_marked[scc_count] = 1;
                    } while (w != cur->q);
                    scc_size[scc_count] = sz;
                    scc_count++;
                }
                if (dfs_ptr > 0) {
                    if (lowlink[cur->q] < lowlink[dfs[dfs_ptr-1].q])
                        lowlink[dfs[dfs_ptr-1].q] = lowlink[cur->q];
                }
            }
        }
    }

    uint16_t livelock_states = 0;
    for (uint16_t s = 0; s < scc_count; s++) {
        if (!scc_has_marked[s] && scc_size[s] > 0) {
            /* A livelock is a cycle: check if SCC has internal transitions */
            int has_internal = 0;
            for (uint16_t q = 0; q < G->nstates && !has_internal; q++) {
                if (scc_id[q] != s) continue;
                for (uint16_t ev = 0; ev < G->nevents && !has_internal; ev++) {
                    uint16_t dst = G->delta[q][ev];
                    if (dst != DES_UNDEF_STATE && scc_id[dst] == s) {
                        /* Internal transition found - cycle exists */
                        has_internal = 1;
                    }
                }
            }
            if (!has_internal) continue; /* Trivial SCC = deadlock, not livelock */

            int can_escape = 0;
            for (uint16_t q = 0; q < G->nstates && !can_escape; q++) {
                if (scc_id[q] != s) continue;
                for (uint16_t ev = 0; ev < G->nevents && !can_escape; ev++) {
                    uint16_t dst = G->delta[q][ev];
                    if (dst != DES_UNDEF_STATE && scc_id[dst] != s) can_escape = 1;
                }
            }
            if (!can_escape) livelock_states += scc_size[s];
        }
    }
    return livelock_states;
}

void des_automaton_destroy(des_automaton_t *G)
{
    if (!G) return;
    memset(G, 0, sizeof(*G));
}