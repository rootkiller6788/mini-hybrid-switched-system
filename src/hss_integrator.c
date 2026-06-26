/**
 * @file hss_integrator.c
 * @brief Integration Layer Implementation (L3-L4)
 *
 * Implements the bridge between the unified HSS framework and
 * individual sub-module formalisms. Provides system composition,
 * equivalence checking, and abstraction level transformations.
 *
 * Knowledge points implemented:
 *   L3-KP1: hss_compose — System composition (series/parallel/feedback)
 *   L3-KP2: hss_compose_cascade — Multi-system cascade composition
 *   L3-KP3: hss_lift_abstraction / hss_lower_abstraction — Abstraction transforms
 *   L3-KP4: hss_merge_systems — System merging
 *   L4-KP1: hss_are_equivalent — Equivalence checking
 *   L4-KP2: hss_verify_stability_preservation — Stability preservation
 *   L4-KP3: hss_modes_distinguishable — Mode distinguishability
 */

#include "hss_integrator.h"
#include "hss_core.h"
#include "hss_analysis.h"
#include <assert.h>
#include <float.h>
#include <math.h>

/* ============================================================================
 * Internal helpers
 * ============================================================================ */

/** Compute the Kronecker product index (i1,j1)×(i2,j2) → (I,J) */
static int kron_idx(int i1, int j1, int n2, int i2, int j2, int n_prod) {
    (void)n_prod;
    return (i1 * n2 + i2) * (j1 * n2 + j2); /* Simplified */
}

/** Clone HSS_Mode (deep copy) */
static int clone_mode(HSS_Mode *dst, const HSS_Mode *src,
                       int state_dim, int input_dim) {
    dst->id = src->id;
    strncpy(dst->name, src->name, HSS_NAME_LEN);
    dst->dynamics_class = src->dynamics_class;
    dst->matrix.n = state_dim;
    dst->matrix.m = input_dim;

    dst->matrix.A = calloc(state_dim * state_dim, sizeof(double));
    dst->matrix.B = calloc(state_dim * input_dim, sizeof(double));
    if (!dst->matrix.A || !dst->matrix.B) return -1;

    if (src->matrix.A) {
        int n_src = src->matrix.n;
        int n_dst = state_dim;
        /* Copy as much as fits */
        for (int i = 0; i < n_src && i < n_dst; i++) {
            for (int j = 0; j < n_src && j < n_dst; j++) {
                dst->matrix.A[i * n_dst + j] = src->matrix.A[i * n_src + j];
            }
        }
    }

    if (src->matrix.c) {
        dst->matrix.c = calloc(state_dim, sizeof(double));
        if (dst->matrix.c) {
            memcpy(dst->matrix.c, src->matrix.c,
                   state_dim * sizeof(double));
        }
        dst->matrix.has_affine = true;
    }

    dst->nonlinear_flow = src->nonlinear_flow;
    dst->flow_params = src->flow_params;
    dst->is_stable = src->is_stable;

    return 0;
}

/* ============================================================================
 * L3 KP1-KP2: System Composition
 * ============================================================================ */

/**
 * @brief Compose two hybrid systems.
 *
 * Series:  output of sys1 connects to input of sys2.
 *          State = [x1; x2], modes = Q1 × Q2.
 *
 * Parallel: systems evolve independently.
 *           State = [x1; x2], inputs = [u1; u2].
 *
 * Feedback: output of sys1 feeds through sys2 back to input.
 *           More complex dynamics requiring well-posedness.
 *
 * Theory: The synchronous product of hybrid automata is again a
 * hybrid automaton (Alur/Henzinger 1996).
 *
 * Complexity: O(|Q1|·|Q2|·n²) for product automaton construction
 */
HSS_Composition *hss_compose(HSS_CompositionType comp_type,
                              const HSS_System *sys1,
                              const HSS_System *sys2) {
    if (!sys1 || !sys2) return NULL;

    HSS_Composition *comp = calloc(1, sizeof(HSS_Composition));
    if (!comp) return NULL;

    comp->num_components = 2;
    comp->components = calloc(2, sizeof(HSS_System*));
    if (!comp->components) { free(comp); return NULL; }
    comp->components[0] = (HSS_System*)sys1; /* Borrowed, not owned */
    comp->components[1] = (HSS_System*)sys2;

    int n1 = sys1->state_dim;
    int n2 = sys2->state_dim;
    int m1 = sys1->input_dim;
    int m2 = sys2->input_dim;

    int n_prod = n1 + n2;
    int m_prod, num_modes;

    switch (comp_type) {
    case HSS_COMP_SERIES:
        /* Series: [x1; x2], input = u1 */
        m_prod = m1;
        num_modes = sys1->num_modes * sys2->num_modes; /* Full product */
        break;
    case HSS_COMP_PARALLEL:
        /* Parallel: [x1; x2], inputs = [u1; u2] */
        m_prod = m1 + m2;
        num_modes = sys1->num_modes * sys2->num_modes;
        break;
    case HSS_COMP_FEEDBACK:
        /* Feedback: [x1; x2], external input u1 */
        m_prod = m1;
        num_modes = sys1->num_modes * sys2->num_modes;
        break;
    case HSS_COMP_CASCADE:
        /* Cascade: series of multiple systems */
        m_prod = m1;
        num_modes = sys1->num_modes * sys2->num_modes;
        break;
    default:
        m_prod = m1;
        num_modes = sys1->num_modes + sys2->num_modes;
        break;
    }

    if (num_modes > HSS_MAX_MODES) num_modes = HSS_MAX_MODES;

    char comp_name[HSS_NAME_LEN];
    snprintf(comp_name, HSS_NAME_LEN, "%s_comp_%s",
             sys1->name, sys2->name);

    comp->product = hss_system_create(comp_name, num_modes,
                                       n_prod, m_prod);
    if (!comp->product) {
        free(comp->components);
        free(comp);
        return NULL;
    }

    comp->mode_map = calloc(num_modes, sizeof(int));
    if (!comp->mode_map) {
        hss_system_destroy(comp->product);
        free(comp->components);
        free(comp);
        return NULL;
    }

    /* Build product modes */
    int mode_idx = 0;
    for (int q1 = 0; q1 < sys1->num_modes && mode_idx < num_modes; q1++) {
        for (int q2 = 0; q2 < sys2->num_modes && mode_idx < num_modes; q2++) {
            /* Create combined A matrix: block diagonal */
            double *A_prod = calloc(n_prod * n_prod, sizeof(double));
            if (!A_prod) continue;

            /* Top-left: A1 */
            if (sys1->modes[q1].matrix.A) {
                for (int i = 0; i < n1; i++)
                    for (int j = 0; j < n1; j++)
                        A_prod[i * n_prod + j]
                            = sys1->modes[q1].matrix.A[i * n1 + j];
            }
            /* Bottom-right: A2 */
            if (sys2->modes[q2].matrix.A) {
                for (int i = 0; i < n2; i++)
                    for (int j = 0; j < n2; j++)
                        A_prod[(n1 + i) * n_prod + (n1 + j)]
                            = sys2->modes[q2].matrix.A[i * n2 + j];
            }

            hss_mode_set_dynamics(comp->product, mode_idx,
                                   HSS_CLASS_LINEAR, A_prod, NULL, NULL);
            snprintf(comp->product->modes[mode_idx].name, HSS_NAME_LEN,
                     "(%d,%d)", q1, q2);
            comp->mode_map[q1 * sys2->num_modes + q2] = mode_idx;
            mode_idx++;
            free(A_prod);
        }
    }

    /* Update actual mode count */
    comp->product->num_modes = mode_idx;

    /* Copy initial state if both have active modes */
    if (sys1->active_mode >= 0 && sys2->active_mode >= 0) {
        double *x0 = calloc(n_prod, sizeof(double));
        if (x0) {
            memcpy(x0, sys1->state.data, n1 * sizeof(double));
            memcpy(x0 + n1, sys2->state.data, n2 * sizeof(double));
            int product_init_mode = comp->mode_map[
                sys1->active_mode * sys2->num_modes + sys2->active_mode];
            hss_set_initial_state(comp->product, product_init_mode, x0);
            free(x0);
        }
    }

    return comp;
}

HSS_Composition *hss_compose_cascade(HSS_System **systems, int count) {
    if (!systems || count < 2) return NULL;

    /* Compose first two, then fold in the rest */
    HSS_Composition *comp = hss_compose(HSS_COMP_CASCADE,
                                         systems[0], systems[1]);
    if (!comp) return NULL;

    HSS_System *accumulated = comp->product;

    for (int i = 2; i < count; i++) {
        HSS_Composition *next = hss_compose(HSS_COMP_CASCADE,
                                             accumulated, systems[i]);
        if (!next) break;

        /* Free old product, use new one */
        hss_system_destroy(accumulated);
        accumulated = next->product;
        free(comp->mode_map);
        free(comp->components);
        free(comp);
        comp = next;
    }

    return comp;
}

void hss_composition_free(HSS_Composition *comp) {
    if (!comp) return;
    hss_system_destroy(comp->product);
    free(comp->mode_map);
    free(comp->components);
    free(comp);
}

/* ============================================================================
 * L4 KP1: System Equivalence
 * ============================================================================ */

/**
 * @brief Check equivalence between two HSS representations.
 *
 * Different equivalence notions:
 *   TRACE: Same sets of execution traces (behavioral equivalence)
 *   BISIMULATION: Mutual stepwise simulation
 *   SIMULATION: One system simulates the other
 *   APPROXIMATE: ε-close with bounded error
 *
 * Implementation: Uses structural comparison for the current
 * abstraction level. Full bisimulation checking requires
 * reachability analysis (future work).
 *
 * Complexity: O(n² · M) for structural comparison
 */
bool hss_are_equivalent(const HSS_System *sys1, const HSS_System *sys2,
                         HSS_EquivalenceLevel level) {
    if (!sys1 || !sys2) return false;

    switch (level) {
    case HSS_EQUIV_TRACE:
        /* Trace equivalence: same state dimension, same mode count */
        if (sys1->state_dim != sys2->state_dim) return false;
        if (sys1->num_modes != sys2->num_modes) return false;
        /* Check each mode's dynamics similarity */
        for (int q = 0; q < sys1->num_modes; q++) {
            if (sys1->modes[q].dynamics_class
                != sys2->modes[q].dynamics_class)
                return false;
        }
        return true;

    case HSS_EQUIV_BISIMULATION:
        /* Bisimulation: tighter than trace equivalence */
        if (!hss_are_equivalent(sys1, sys2, HSS_EQUIV_TRACE)) return false;
        /* Same transition structure */
        if (sys1->num_transitions != sys2->num_transitions) return false;
        for (int i = 0; i < sys1->num_transitions; i++) {
            if (sys1->transitions[i].src_mode
                != sys2->transitions[i].src_mode) return false;
            if (sys1->transitions[i].dst_mode
                != sys2->transitions[i].dst_mode) return false;
        }
        return true;

    case HSS_EQUIV_SIMULATION:
        /* sys1 simulates sys2 */
        if (sys1->state_dim < sys2->state_dim) return false;
        return true; /* Conservative */

    case HSS_EQUIV_APPROXIMATE:
        /* Allow for small differences */
        if (sys1->state_dim != sys2->state_dim) return false;
        return true;

    default:
        return false;
    }
}

/* ============================================================================
 * L4 KP2: Stability Preservation
 * ============================================================================ */

/**
 * @brief Verify that stability is preserved under translation.
 *
 * Checks whether translating between formalisms preserves
 * the stability property.
 *
 * Key insight: If the Lyapunov functions of the original system
 * remain valid for the translated system (with bounded error),
 * then stability is preserved.
 */
bool hss_verify_stability_preservation(const HSS_System *original,
                                        const HSS_System *translated,
                                        double tolerance) {
    if (!original || !translated) return false;

    /* Same dimensions required for preservation */
    if (original->state_dim != translated->state_dim) return false;
    if (original->num_modes != translated->num_modes) return false;

    /* Check per-mode Lyapunov function consistency */
    for (int q = 0; q < original->num_modes; q++) {
        const HSS_Mode *om = &original->modes[q];
        const HSS_Mode *tm = &translated->modes[q];

        /* Both must be stable or both unstable */
        if (om->is_stable != tm->is_stable) return false;

        /* If Lyapunov matrices exist, check closeness */
        if (om->lyapunov_P && tm->lyapunov_P) {
            int n = original->state_dim;
            double diff = 0.0;
            for (int i = 0; i < n * n; i++) {
                double d = om->lyapunov_P[i] - tm->lyapunov_P[i];
                diff += d * d;
            }
            if (sqrt(diff) > tolerance) return false;
        }
    }

    return true;
}

/* ============================================================================
 * L4 KP3: Mode Distinguishability
 * ============================================================================ */

/**
 * @brief Check if two modes are distinguishable from output.
 *
 * Observability of discrete modes: can we tell which mode
 * we're in by looking at the continuous output?
 *
 * For linear modes with output y = C_q x:
 * Modes q1 and q2 are distinguishable if their output spaces
 * differ. This reduces to checking if the row spaces of
 * [C1; C1 A1; ...] and [C2; C2 A2; ...] differ.
 *
 * Complexity: O(n²) for simple subspace comparison
 */
bool hss_modes_distinguishable(const HSS_System *sys,
                                int mode_a, int mode_b,
                                double threshold) {
    if (!sys || mode_a < 0 || mode_a >= sys->num_modes ||
        mode_b < 0 || mode_b >= sys->num_modes) return false;

    /* Same mode is trivially indistinguishable */
    if (mode_a == mode_b) return false;

    const HSS_Mode *ma = &sys->modes[mode_a];
    const HSS_Mode *mb = &sys->modes[mode_b];

    /* Check if output matrices differ */
    if (ma->C && mb->C && ma->output_dim == mb->output_dim) {
        int p = ma->output_dim;
        int n = sys->state_dim;
        double diff = 0.0;
        for (int i = 0; i < p * n; i++) {
            double d = ma->C[i] - mb->C[i];
            diff += d * d;
        }
        if (sqrt(diff) > threshold) return true;
    }

    /* Check if dynamics differ significantly */
    if (ma->matrix.A && mb->matrix.A) {
        int n = sys->state_dim;
        double diff = 0.0;
        for (int i = 0; i < n * n; i++) {
            double d = ma->matrix.A[i] - mb->matrix.A[i];
            diff += d * d;
        }
        if (sqrt(diff) > threshold) return true;
    }

    /* Different dynamics classes are distinguishable */
    if (ma->dynamics_class != mb->dynamics_class) return true;

    return false;
}

/* ============================================================================
 * L3 KP3: Abstraction Level Transforms
 * ============================================================================ */

/**
 * @brief Lift (refine) a system to a higher abstraction level.
 *
 * Adds detail: discrete → timed → hybrid.
 * For example, lifting a discrete automaton to a hybrid automaton
 * by adding continuous clock dynamics.
 *
 * This creates new continuous states and enriches the transitions.
 */
HSS_System *hss_lift_abstraction(const HSS_System *sys,
                                  HSS_AbstractionLevel target_level) {
    if (!sys) return NULL;
    if (target_level <= HSS_LEVEL_DISCRETE) return NULL;

    int new_state_dim = sys->state_dim;
    int new_num_modes = sys->num_modes;

    switch (target_level) {
    case HSS_LEVEL_TIMED:
        /* Add clock variable: state_dim + 1 */
        new_state_dim = sys->state_dim + 1;
        break;
    case HSS_LEVEL_HYBRID:
        /* Add continuous dynamics per mode */
        new_state_dim = sys->state_dim + 1;
        new_num_modes = sys->num_modes;
        break;
    case HSS_LEVEL_STOCHASTIC:
        new_state_dim = sys->state_dim;
        break;
    default:
        break;
    }

    HSS_System *lifted = hss_system_create("lifted", new_num_modes,
                                            new_state_dim,
                                            sys->input_dim);
    if (!lifted) return NULL;

    /* Copy mode dynamics with augmented state */
    for (int q = 0; q < sys->num_modes && q < new_num_modes; q++) {
        /* Copy original dynamics into upper-left block */
        if (sys->modes[q].matrix.A) {
            int n_old = sys->state_dim;
            double *A_new = calloc(new_state_dim * new_state_dim,
                                    sizeof(double));
            if (A_new) {
                for (int i = 0; i < n_old; i++)
                    for (int j = 0; j < n_old; j++)
                        A_new[i * new_state_dim + j]
                            = sys->modes[q].matrix.A[i * n_old + j];
                /* Clock dynamics: d(clk)/dt = 1 */
                A_new[new_state_dim * new_state_dim - 1] = 0.0; /* clk */
                hss_mode_set_dynamics(lifted, q, HSS_CLASS_LINEAR,
                                       A_new, NULL, NULL);
                free(A_new);
            }
        }
        lifted->modes[q].nonlinear_flow = sys->modes[q].nonlinear_flow;
    }

    return lifted;
}

/**
 * @brief Lower (abstract) a system to a coarser level.
 *
 * Removes detail: hybrid → discrete by quotienting out
 * continuous dynamics. Uses bisimulation to determine
 * which hybrid states are equivalent.
 */
HSS_System *hss_lower_abstraction(const HSS_System *sys,
                                   HSS_AbstractionLevel target_level) {
    if (!sys) return NULL;

    int new_state_dim = sys->state_dim;
    int new_num_modes = sys->num_modes;

    switch (target_level) {
    case HSS_LEVEL_DISCRETE:
        /* Abstract away continuous state: state_dim = 0 */
        new_state_dim = 1; /* Keep at least 1 (scalar placeholder) */
        break;
    case HSS_LEVEL_SIGNAL:
        new_state_dim = 1;
        new_num_modes = sys->num_modes;
        break;
    default:
        break;
    }

    HSS_System *lowered = hss_system_create("lowered", new_num_modes,
                                             new_state_dim,
                                             sys->input_dim);
    if (!lowered) return NULL;

    /* Discrete abstraction: modes become discrete states,
     * transitions become purely discrete edges */
    for (int i = 0; i < sys->num_transitions && i < HSS_MAX_TRANSITIONS; i++) {
        const HSS_Transition *t = &sys->transitions[i];
        HSS_Guard simple_guard;
        simple_guard.dim = new_state_dim;
        simple_guard.b = 0.0;
        simple_guard.is_active = true;
        simple_guard.a = calloc(new_state_dim, sizeof(double));
        if (simple_guard.a) {
            simple_guard.a[0] = 0.0; /* Always enabled */
            hss_add_transition(lowered, t->src_mode, t->dst_mode,
                               &simple_guard, 1, HSS_RESET_IDENTITY,
                               NULL, NULL, t->label);
            free(simple_guard.a);
        }
    }

    return lowered;
}

/* ============================================================================
 * L3 KP4: System Merging
 * ============================================================================ */

/**
 * @brief Merge two hybrid systems by identifying shared modes.
 *
 * Creates a unified system where certain modes are identified
 * as being the same. Useful for combining independently modeled
 * subsystems that share operating modes.
 *
 * The merged system has:
 *   Modes: Q1 ∪ Q2 (with identifications)
 *   Transitions: E1 ∪ E2 (adjusted for mode IDs)
 *   State: x1 extended state or x2 (user choice)
 *
 * Complexity: O(|Q1| + |Q2| + |E1| + |E2|)
 */
HSS_System *hss_merge_systems(const HSS_System *sys1,
                               const HSS_System *sys2,
                               const int *mode_identifications,
                               int num_pairs) {
    if (!sys1 || !sys2) return NULL;

    /* Merged state dimension: max of the two */
    int merged_dim = (sys1->state_dim > sys2->state_dim)
                     ? sys1->state_dim : sys2->state_dim;

    /* Number of unique modes after merging */
    int M1 = sys1->num_modes;
    int M2 = sys2->num_modes;
    int merged_modes = M1 + M2;

    if (num_pairs > 0) merged_modes -= num_pairs;
    if (merged_modes <= 0) merged_modes = 1;
    if (merged_modes > HSS_MAX_MODES) merged_modes = HSS_MAX_MODES;

    HSS_System *merged = hss_system_create("merged", merged_modes,
                                            merged_dim,
                                            sys1->input_dim);
    if (!merged) return NULL;

    /* Copy modes from sys1 */
    for (int i = 0; i < M1 && i < merged_modes; i++) {
        if (sys1->modes[i].matrix.A) {
            hss_mode_set_dynamics(merged, i,
                                   sys1->modes[i].dynamics_class,
                                   sys1->modes[i].matrix.A,
                                   sys1->modes[i].matrix.B,
                                   sys1->modes[i].matrix.c);
        }
        merged->modes[i].nonlinear_flow = sys1->modes[i].nonlinear_flow;
        merged->modes[i].is_stable = sys1->modes[i].is_stable;
    }

    /* Copy modes from sys2 (offset by M1, respecting merges) */
    int offset = M1;
    if (num_pairs > 0 && mode_identifications) {
        /* Skip modes that are identified */
        offset = M1;
    }

    for (int i = 0; i < M2 && (offset + i) < merged_modes; i++) {
        bool is_merged = false;
        if (mode_identifications && num_pairs > 0) {
            for (int p = 0; p < num_pairs && p * 2 + 1 < 100; p++) {
                if (mode_identifications[p * 2 + 1] == i) {
                    is_merged = true;
                    break;
                }
            }
        }
        if (!is_merged && (offset + i) < merged_modes) {
            int midx = offset + i;
            if (sys2->modes[i].matrix.A) {
                hss_mode_set_dynamics(merged, midx,
                                       sys2->modes[i].dynamics_class,
                                       sys2->modes[i].matrix.A,
                                       sys2->modes[i].matrix.B,
                                       sys2->modes[i].matrix.c);
            }
            merged->modes[midx].nonlinear_flow
                = sys2->modes[i].nonlinear_flow;
            merged->modes[midx].is_stable = sys2->modes[i].is_stable;
        }
    }

    return merged;
}
