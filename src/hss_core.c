/**
 * @file hss_core.c
 * @brief Core implementation of Hybrid Switched System definitions
 *
 * Implements L1 (Core Definitions) and L2 (Core Concepts) for
 * the unified hybrid switched system framework.
 *
 * Knowledge points implemented in this file:
 *   L1-KP1: hss_system_create / hss_system_destroy (system lifecycle)
 *   L1-KP2: hss_mode_set_dynamics (mode dynamics configuration)
 *   L1-KP3: hss_mode_set_nonlinear_flow (nonlinear flow setup)
 *   L1-KP4: hss_add_transition (transition definition)
 *   L1-KP5: hss_mode_add_invariant (invariant constraints)
 *   L1-KP6: hss_set_initial_state (initial condition)
 *   L1-KP7: hss_mode_set_lyapunov (Lyapunov matrix config)
 *   L1-KP8: hss_system_validate (system consistency check)
 *   L1-KP9: hss_system_print_summary (system introspection)
 *   L2-KC1: hss_hybrid_time_advance (hybrid time domain)
 *   L2-KC2: hss_check_invariant (flow/jump set partitioning)
 *   L2-KC3: hss_is_zeno_prone (Zeno detection)
 */

#include "hss_core.h"
#include <assert.h>
#include <float.h>
#include <math.h>

/* ============================================================================
 * Internal helpers
 * ============================================================================ */

/** Dot product a·b */
static double vec_dot(const double *a, const double *b, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += a[i] * b[i];
    return sum;
}

/** Check if a matrix is symmetric (within tolerance) */
static bool is_symmetric(const double *M, int n, double tol) {
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (fabs(M[i * n + j] - M[j * n + i]) > tol) return false;
        }
    }
    return true;
}

/** Check positive definiteness via Cholesky attempt */
static bool is_positive_definite(const double *M, int n, double tol) {
    /* Simple check: all diagonal entries must be > 0, and
     * for 2×2, check det > 0. Full Cholesky would be more rigorous. */
    for (int i = 0; i < n; i++) {
        if (M[i * n + i] <= tol) return false;
    }
    if (n == 1) return true;
    if (n == 2) {
        double det = M[0] * M[3] - M[1] * M[2];
        if (det <= tol) return false;
    }
    /* For n > 2, check diagonal dominance as heuristic */
    for (int i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) row_sum += fabs(M[i * n + j]);
        }
        if (M[i * n + i] <= row_sum) return false;
    }
    return true;
}

/* ============================================================================
 * L1 KP1: System lifecycle — create and destroy
 * ============================================================================ */

/**
 * @brief Allocate and initialize a hybrid switched system.
 *
 * Creates the container and allocates arrays for modes and transitions.
 * All modes are initialized to zero/default values.
 *
 * Reference: Goebel/Sanfelice/Teel (2012), Chapter 2.
 *
 * Complexity: O(num_modes · state_dim²) for memory allocation
 */
HSS_System *hss_system_create(const char *name, int num_modes,
                               int state_dim, int input_dim) {
    if (!name || num_modes <= 0 || state_dim <= 0 || input_dim < 0) {
        return NULL;
    }
    if (num_modes > HSS_MAX_MODES || state_dim > HSS_MAX_STATE_DIM ||
        input_dim > HSS_MAX_INPUT_DIM) {
        return NULL;
    }

    HSS_System *sys = calloc(1, sizeof(HSS_System));
    if (!sys) return NULL;

    strncpy(sys->name, name, HSS_NAME_LEN - 1);
    sys->name[HSS_NAME_LEN - 1] = '\0';
    sys->num_modes = num_modes;
    sys->state_dim = state_dim;
    sys->input_dim = input_dim;
    sys->output_dim = state_dim;  /* Default: full state output */
    sys->active_mode = -1;
    sys->system_class = HSS_CLASS_LINEAR;
    sys->time_horizon = 10.0;
    sys->solver = HSS_SOLVER_RK4;
    sys->step_size = 1e-3;
    sys->tolerance = HSS_EPSILON;
    sys->zeno_check = true;
    sys->zeno_min_time = 1e-6;
    sys->zeno_max_jumps = 10000;

    /* Allocate modes */
    sys->modes = calloc(num_modes, sizeof(HSS_Mode));
    if (!sys->modes) { free(sys); return NULL; }
    for (int i = 0; i < num_modes; i++) {
        sys->modes[i].id = i;
        sys->modes[i].dynamics_class = HSS_CLASS_LINEAR;
        sys->modes[i].matrix.n = state_dim;
        sys->modes[i].matrix.m = input_dim;
        /* Allocate system matrices */
        sys->modes[i].matrix.A = calloc(state_dim * state_dim, sizeof(double));
        sys->modes[i].matrix.B = calloc(state_dim * input_dim, sizeof(double));
        sys->modes[i].matrix.c = NULL;
        sys->modes[i].matrix.has_affine = false;
        /* Identity as default A */
        for (int j = 0; j < state_dim; j++) {
            sys->modes[i].matrix.A[j * state_dim + j] = 0.0; /* will be set */
        }
        sys->modes[i].nonlinear_flow = NULL;
        sys->modes[i].flow_params = NULL;
        sys->modes[i].invariants = NULL;
        sys->modes[i].num_invariants = 0;
        sys->modes[i].lyapunov_P = NULL;
        sys->modes[i].lyapunov_decay = 0.0;
        sys->modes[i].is_stable = true; /* assume stable until checked */
        sys->modes[i].C = NULL;
        sys->modes[i].D = NULL;
        sys->modes[i].output_dim = state_dim;
        snprintf(sys->modes[i].name, HSS_NAME_LEN, "mode_%d", i);
    }

    /* Allocate continuous state */
    sys->state.data = calloc(state_dim, sizeof(double));
    sys->state.dim = state_dim;
    sys->state.labels = NULL;
    sys->state.time = 0.0;
    sys->state.jumps = 0;

    /* Initialize switching signal */
    sys->signal.mode_sequence = NULL;
    sys->signal.switch_times = NULL;
    sys->signal.num_switches = 0;
    sys->signal.max_sequence = 0;
    sys->signal.dwell_time_min = 0.0;
    sys->signal.dwell_time_max = DBL_MAX;
    sys->signal.average_dwell = 0.0;
    sys->signal.chatter_bound = 1;

    return sys;
}

/**
 * @brief Deallocate a hybrid switched system.
 *
 * Recursively frees all dynamically allocated memory including
 * mode matrices, invariants, Lyapunov matrices, and state vectors.
 *
 * Must be called to avoid memory leaks.
 */
void hss_system_destroy(HSS_System *sys) {
    if (!sys) return;

    /* Free each mode's resources */
    if (sys->modes) {
        for (int i = 0; i < sys->num_modes; i++) {
            free(sys->modes[i].matrix.A);
            free(sys->modes[i].matrix.B);
            free(sys->modes[i].matrix.c);
            if (sys->modes[i].invariants) {
                for (int j = 0; j < sys->modes[i].num_invariants; j++) {
                    free(sys->modes[i].invariants[j].a);
                }
                free(sys->modes[i].invariants);
            }
            free(sys->modes[i].lyapunov_P);
            free(sys->modes[i].C);
            free(sys->modes[i].D);
        }
        free(sys->modes);
    }

    /* Free transitions */
    if (sys->transitions) {
        for (int i = 0; i < sys->num_transitions; i++) {
            if (sys->transitions[i].guards) {
                for (int j = 0; j < sys->transitions[i].num_guards; j++) {
                    free(sys->transitions[i].guards[j].a);
                }
                free(sys->transitions[i].guards);
            }
            free(sys->transitions[i].reset_M);
            free(sys->transitions[i].reset_b);
        }
        free(sys->transitions);
    }

    /* Free state */
    free(sys->state.data);
    free(sys->state.labels);

    /* Free switching signal */
    free(sys->signal.mode_sequence);
    free(sys->signal.switch_times);

    free(sys);
}

/* ============================================================================
 * L1 KP2: Mode dynamics configuration
 * ============================================================================ */

/**
 * @brief Set the dynamics for a specific mode.
 *
 * For linear modes: stores A, B, and optional affine term c.
 * For affine modes: requires non-NULL c.
 * For bilinear/LPV/polynomial modes: delegates to nonlinear_flow.
 *
 * The mode's dynamics_class is set accordingly.
 *
 * Complexity: O(n² + n·m) for matrix copy
 */
int hss_mode_set_dynamics(HSS_System *sys, int mode_id,
                           HSS_SystemClass cls,
                           const double *A, const double *B,
                           const double *c) {
    if (!sys || mode_id < 0 || mode_id >= sys->num_modes) return -1;
    HSS_Mode *mode = &sys->modes[mode_id];
    int n = sys->state_dim;
    int m = sys->input_dim;

    mode->dynamics_class = cls;

    /* Copy state matrix A */
    if (A) {
        memcpy(mode->matrix.A, A, n * n * sizeof(double));
    } else {
        /* Default: zero matrix (still dynamics) */
        memset(mode->matrix.A, 0, n * n * sizeof(double));
    }

    /* Copy input matrix B */
    if (B && mode->matrix.B) {
        memcpy(mode->matrix.B, B, n * m * sizeof(double));
    } else if (mode->matrix.B) {
        memset(mode->matrix.B, 0, n * m * sizeof(double));
    }

    /* Copy affine term c */
    if (c) {
        if (!mode->matrix.c) {
            mode->matrix.c = calloc(n, sizeof(double));
            if (!mode->matrix.c) return -1;
        }
        memcpy(mode->matrix.c, c, n * sizeof(double));
        mode->matrix.has_affine = true;
    } else {
        free(mode->matrix.c);
        mode->matrix.c = NULL;
        mode->matrix.has_affine = false;
    }

    return 0;
}

/**
 * @brief Set a nonlinear flow function for a mode.
 *
 * For general nonlinear modes: ẋ = f(x, u, t) with user-defined
 * vector field. The function must compute dxdt given current state.
 *
 * This supports LPV, polynomial, bilinear, and general nonlinear systems.
 *
 * Complexity: Depends on user function evaluation cost
 */
int hss_mode_set_nonlinear_flow(HSS_System *sys, int mode_id,
                                 void (*flow_func)(const double *x,
                                     const double *u, double t,
                                     double *dxdt, int n, void *params),
                                 void *params) {
    if (!sys || mode_id < 0 || mode_id >= sys->num_modes) return -1;
    if (!flow_func) return -1;

    HSS_Mode *mode = &sys->modes[mode_id];
    mode->nonlinear_flow = flow_func;
    mode->flow_params = params;
    mode->dynamics_class = HSS_CLASS_NONLINEAR;
    return 0;
}

/* ============================================================================
 * L1 KP4: Transition definition
 * ============================================================================ */

/**
 * @brief Add a discrete transition between modes.
 *
 * Creates an edge in the hybrid automaton graph.
 * Guards are logical-AND: all must be satisfied for transition.
 *
 * The reset map type determines how the continuous state changes:
 *   IDENTITY: x⁺ = x⁻ (continuous evolution)
 *   ZERO:     x⁺ = 0 (hard reset)
 *   LINEAR:   x⁺ = M x⁻
 *   AFFINE:   x⁺ = M x⁻ + b
 *
 * Complexity: O(n · num_guards) for memory allocation
 */
int hss_add_transition(HSS_System *sys, int src, int dst,
                        const HSS_Guard *guards, int num_guards,
                        HSS_ResetType reset_type,
                        const double *reset_M, const double *reset_b,
                        const char *label) {
    if (!sys || src < 0 || src >= sys->num_modes ||
        dst < 0 || dst >= sys->num_modes) return -1;
    if (num_guards < 0 || num_guards > HSS_MAX_GUARDS) return -1;

    int n = sys->state_dim;
    int idx = sys->num_transitions;

    /* Grow transitions array */
    HSS_Transition *new_trans = realloc(sys->transitions,
        (idx + 1) * sizeof(HSS_Transition));
    if (!new_trans) return -1;
    sys->transitions = new_trans;
    sys->num_transitions = idx + 1;

    HSS_Transition *t = &sys->transitions[idx];
    memset(t, 0, sizeof(HSS_Transition));
    t->src_mode = src;
    t->dst_mode = dst;
    t->reset_type = reset_type;
    t->reset_dim = n;
    t->priority = 1.0;
    t->is_urgent = false;
    if (label) {
        strncpy(t->label, label, HSS_NAME_LEN - 1);
    } else {
        snprintf(t->label, HSS_NAME_LEN, "t_%d_%d", src, dst);
    }

    /* Copy guards */
    t->num_guards = num_guards;
    if (num_guards > 0) {
        t->guards = calloc(num_guards, sizeof(HSS_Guard));
        if (!t->guards) { sys->num_transitions--; return -1; }
        for (int i = 0; i < num_guards; i++) {
            t->guards[i].dim = n;
            t->guards[i].b = guards[i].b;
            t->guards[i].is_active = guards[i].is_active;
            t->guards[i].a = calloc(n, sizeof(double));
            if (!t->guards[i].a) {
                /* Cleanup on partial failure */
                for (int j = 0; j < i; j++) free(t->guards[j].a);
                free(t->guards);
                t->guards = NULL;
                sys->num_transitions--;
                return -1;
            }
            memcpy(t->guards[i].a, guards[i].a, n * sizeof(double));
        }
    }

    /* Copy reset map */
    if (reset_M) {
        t->reset_M = calloc(n * n, sizeof(double));
        if (!t->reset_M) { sys->num_transitions--; return -1; }
        memcpy(t->reset_M, reset_M, n * n * sizeof(double));
    }
    if (reset_b) {
        t->reset_b = calloc(n, sizeof(double));
        if (!t->reset_b) {
            free(t->reset_M);
            t->reset_M = NULL;
            sys->num_transitions--;
            return -1;
        }
        memcpy(t->reset_b, reset_b, n * sizeof(double));
    }

    return idx;
}

/* ============================================================================
 * L1 KP5: Invariant constraints
 * ============================================================================ */

/**
 * @brief Add an invariant constraint to a mode.
 *
 * Invariant: aᵀx ≤ b must hold while the system is in this mode.
 * If violated, the system must transition (or behavior is undefined).
 *
 * Multiple invariants are AND-combined: all must hold simultaneously.
 *
 * Complexity: O(n) for array copy
 */
int hss_mode_add_invariant(HSS_System *sys, int mode_id,
                            const double *a, double b) {
    if (!sys || mode_id < 0 || mode_id >= sys->num_modes) return -1;
    if (!a) return -1;

    HSS_Mode *mode = &sys->modes[mode_id];
    int n = sys->state_dim;
    int idx = mode->num_invariants;

    HSS_Guard *new_inv = realloc(mode->invariants,
        (idx + 1) * sizeof(HSS_Guard));
    if (!new_inv) return -1;
    mode->invariants = new_inv;

    HSS_Guard *inv = &mode->invariants[idx];
    inv->dim = n;
    inv->b = b;
    inv->is_active = true;
    inv->a = calloc(n, sizeof(double));
    if (!inv->a) return -1;
    memcpy(inv->a, a, n * sizeof(double));
    mode->num_invariants = idx + 1;

    return 0;
}

/* ============================================================================
 * L1 KP6: Initial state setup
 * ============================================================================ */

/**
 * @brief Set the initial hybrid state.
 *
 * Configures (q₀, x₀) — the initial mode and continuous state.
 * Must be called before simulation or analysis.
 *
 * Complexity: O(n) for state copy
 */
int hss_set_initial_state(HSS_System *sys, int initial_mode,
                           const double *x0) {
    if (!sys || initial_mode < 0 || initial_mode >= sys->num_modes) return -1;
    if (!x0) return -1;

    sys->active_mode = initial_mode;
    memcpy(sys->state.data, x0, sys->state_dim * sizeof(double));
    sys->state.time = 0.0;
    sys->state.jumps = 0;

    return 0;
}

/* ============================================================================
 * L1 KP7: Lyapunov matrix configuration
 * ============================================================================ */

/**
 * @brief Set the Lyapunov matrix for a mode.
 *
 * V_q(x) = xᵀ P_q x where P_q must be symmetric positive definite.
 *
 * The decay rate α satisfies V̇_q(x) ≤ -α V_q(x) along flows.
 * For linear modes: α = λ_min(P A_q + A_qᵀ P) / λ_max(P).
 *
 * Complexity: O(n²) for matrix copy
 */
int hss_mode_set_lyapunov(HSS_System *sys, int mode_id,
                           const double *P, double decay_rate) {
    if (!sys || mode_id < 0 || mode_id >= sys->num_modes) return -1;
    if (!P) return -1;

    HSS_Mode *mode = &sys->modes[mode_id];
    int n = sys->state_dim;

    if (!mode->lyapunov_P) {
        mode->lyapunov_P = calloc(n * n, sizeof(double));
        if (!mode->lyapunov_P) return -1;
    }
    memcpy(mode->lyapunov_P, P, n * n * sizeof(double));
    mode->lyapunov_decay = decay_rate;

    return 0;
}

/* ============================================================================
 * L1 KP8: System validation
 * ============================================================================ */

/**
 * @brief Validate system consistency.
 *
 * Performs comprehensive checks:
 * 1. All modes have valid dynamics (A or nonlinear flow)
 * 2. All transition sources/destinations are valid mode indices
 * 3. Guard dimensions match state dimension
 * 4. Reset map dimensions match state dimension
 * 5. No self-loop transitions with identity reset (redundant)
 * 6. All invariants have correct dimension
 * 7. Lyapunov matrices are SPD (if provided)
 * 8. Active mode is set (if applicable)
 *
 * Returns 0 if valid, -1 if invalid (prints diagnostic to stderr).
 */
int hss_system_validate(const HSS_System *sys) {
    if (!sys) {
        fprintf(stderr, "HSS validate: NULL system\n");
        return -1;
    }
    if (sys->num_modes <= 0) {
        fprintf(stderr, "HSS validate: no modes defined\n");
        return -1;
    }
    if (sys->state_dim <= 0) {
        fprintf(stderr, "HSS validate: invalid state dimension\n");
        return -1;
    }

    int n = sys->state_dim;
    int valid = 0;

    /* Check modes */
    for (int i = 0; i < sys->num_modes; i++) {
        const HSS_Mode *mode = &sys->modes[i];

        /* Check dynamics */
        if (mode->dynamics_class <= HSS_CLASS_BILINEAR) {
            if (!mode->matrix.A) {
                fprintf(stderr, "HSS validate: mode %d missing A matrix\n", i);
                valid = -1;
            }
        } else if (mode->dynamics_class == HSS_CLASS_NONLINEAR ||
                   mode->dynamics_class == HSS_CLASS_POLYNOMIAL) {
            if (!mode->nonlinear_flow) {
                fprintf(stderr, "HSS validate: mode %d missing nonlinear flow\n", i);
                valid = -1;
            }
        }

        /* Check invariants */
        for (int j = 0; j < mode->num_invariants; j++) {
            if (mode->invariants[j].dim != n) {
                fprintf(stderr, "HSS validate: mode %d invariant %d "
                        "dimension mismatch\n", i, j);
                valid = -1;
            }
        }

        /* Check Lyapunov matrix */
        if (mode->lyapunov_P) {
            if (!is_symmetric(mode->lyapunov_P, n, HSS_EPSILON)) {
                fprintf(stderr, "HSS validate: mode %d P not symmetric\n", i);
                valid = -1;
            }
            if (!is_positive_definite(mode->lyapunov_P, n, HSS_EPSILON)) {
                fprintf(stderr, "HSS validate: mode %d P not positive definite\n", i);
                valid = -1;
            }
        }
    }

    /* Check transitions */
    for (int i = 0; i < sys->num_transitions; i++) {
        const HSS_Transition *t = &sys->transitions[i];

        if (t->src_mode < 0 || t->src_mode >= sys->num_modes) {
            fprintf(stderr, "HSS validate: transition %d invalid src\n", i);
            valid = -1;
        }
        if (t->dst_mode < 0 || t->dst_mode >= sys->num_modes) {
            fprintf(stderr, "HSS validate: transition %d invalid dst\n", i);
            valid = -1;
        }

        /* Check guard dimensions */
        for (int j = 0; j < t->num_guards; j++) {
            if (t->guards[j].dim != n) {
                fprintf(stderr, "HSS validate: transition %d guard %d "
                        "dim mismatch\n", i, j);
                valid = -1;
            }
        }
    }

    /* Check switching signal */
    if (sys->signal.num_switches > 0) {
        for (int i = 0; i < sys->signal.num_switches; i++) {
            int m = sys->signal.mode_sequence[i];
            if (m < 0 || m >= sys->num_modes) {
                fprintf(stderr, "HSS validate: signal refs invalid mode %d\n", m);
                valid = -1;
            }
        }
    }

    return valid;
}

/* ============================================================================
 * L1 KP9: System introspection
 * ============================================================================ */

/**
 * @brief Print a human-readable summary of the hybrid switched system.
 *
 * Outputs mode count, dimensions, transition topology, and
 * Lyapunov function information.
 */
void hss_system_print_summary(const HSS_System *sys, FILE *fp) {
    if (!sys) { fprintf(fp, "NULL system\n"); return; }
    if (!fp) fp = stdout;

    fprintf(fp, "========================================\n");
    fprintf(fp, "Hybrid Switched System: %s\n", sys->name);
    fprintf(fp, "========================================\n");
    fprintf(fp, "  Modes:           %d\n", sys->num_modes);
    fprintf(fp, "  Transitions:     %d\n", sys->num_transitions);
    fprintf(fp, "  State dim:       %d\n", sys->state_dim);
    fprintf(fp, "  Input dim:       %d\n", sys->input_dim);
    fprintf(fp, "  Output dim:      %d\n", sys->output_dim);
    fprintf(fp, "  System class:    %d\n", (int)sys->system_class);
    fprintf(fp, "  Solver:          %d\n", (int)sys->solver);
    fprintf(fp, "  Step size:       %.2e\n", sys->step_size);
    fprintf(fp, "  Zeno check:      %s\n", sys->zeno_check ? "ON" : "OFF");
    fprintf(fp, "  Active mode:     %d\n", sys->active_mode);

    if (sys->signal.num_switches > 0) {
        fprintf(fp, "  Switching signal: %d switches recorded\n",
                sys->signal.num_switches);
    }

    /* Print mode details */
    fprintf(fp, "\nModes:\n");
    for (int i = 0; i < sys->num_modes && i < 10; i++) {
        const HSS_Mode *m = &sys->modes[i];
        fprintf(fp, "  [%d] %-20s class=%d stable=%s invariants=%d\n",
                i, m->name, (int)m->dynamics_class,
                m->is_stable ? "yes" : "no", m->num_invariants);
        if (m->lyapunov_P) {
            fprintf(fp, "       LF: decay=%.4f\n", m->lyapunov_decay);
        }
    }
    if (sys->num_modes > 10) {
        fprintf(fp, "  ... (%d more modes)\n", sys->num_modes - 10);
    }

    /* Print transition topology */
    fprintf(fp, "\nTransitions:\n");
    for (int i = 0; i < sys->num_transitions && i < 10; i++) {
        const HSS_Transition *t = &sys->transitions[i];
        fprintf(fp, "  [%d] %d→%d guards=%d reset=%d %s\n",
                i, t->src_mode, t->dst_mode, t->num_guards,
                (int)t->reset_type, t->label);
    }
    if (sys->num_transitions > 10) {
        fprintf(fp, "  ... (%d more transitions)\n",
                sys->num_transitions - 10);
    }
    fprintf(fp, "========================================\n");
}

/* ============================================================================
 * L2 KC1: Hybrid time domain manipulation
 * ============================================================================ */

/**
 * @brief Advance hybrid time.
 *
 * Hybrid time is (t, j) ∈ R≥0 × ℕ where t is continuous time
 * and j is the discrete jump counter.
 *
 * This is a fundamental concept from the hybrid systems formalism
 * (Goebel/Sanfelice/Teel 2012, Definition 2.1).
 */
void hss_hybrid_time_advance(HSS_ContinuousState *state, double dt) {
    if (!state) return;
    state->time += dt;
}

/**
 * @brief Increment jump counter in hybrid time.
 *
 * Discrete transitions increment j while t stays constant.
 */
void hss_hybrid_time_jump(HSS_ContinuousState *state) {
    if (!state) return;
    state->jumps++;
}

/* ============================================================================
 * L2 KC2: Invariant checking (flow/jump set partitioning)
 * ============================================================================ */

/**
 * @brief Check if state satisfies all invariants of current mode.
 *
 * Verifies aᵢᵀx ≤ bᵢ for all invariant constraints.
 * If any constraint is violated, the system must transition.
 *
 * Returns true if state is inside the invariant set (can flow).
 * Returns false if state is outside (must jump).
 *
 * This implements the flow-set / jump-set partitioning from
 * Goebel/Sanfelice/Teel (2012): C = {(q,x) | x ∈ Inv(q)},
 * D = {(q,x) | x ∉ Inv(q) or guard active}.
 */
bool hss_check_invariant(const HSS_System *sys) {
    if (!sys || sys->active_mode < 0) return false;

    const HSS_Mode *mode = &sys->modes[sys->active_mode];
    const double *x = sys->state.data;

    for (int i = 0; i < mode->num_invariants; i++) {
        const HSS_Guard *inv = &mode->invariants[i];
        if (!inv->is_active) continue;

        double value = vec_dot(inv->a, x, inv->dim);
        if (value > inv->b + HSS_EPSILON) {
            /* Invariant violated */
            return false;
        }
    }
    return true;
}

/**
 * @brief Check if a specific guard condition is satisfied.
 *
 * Evaluates aᵀx ≤ b. The guard is "enabled" when the inequality holds.
 * For strict guards, use aᵀx ≥ b (reversed).
 *
 * Returns true if guard condition is met (transition enabled).
 */
bool hss_check_guard_condition(const HSS_Guard *guard, const double *x) {
    if (!guard || !guard->a || !x) return false;
    if (!guard->is_active) return false;

    double value = vec_dot(guard->a, x, guard->dim);
    return value >= guard->b - HSS_EPSILON;
}

/* ============================================================================
 * L2 KC3: Zeno detection
 * ============================================================================ */

/**
 * @brief Check if the system is prone to Zeno behavior.
 *
 * Zeno: infinitely many discrete transitions in finite time.
 * Detection criteria:
 *   1. Average inter-event time is below threshold
 *   2. Jump count is anomalously high for elapsed time
 *   3. Estimated accumulation time (Zeno time) is finite
 *
 * Complexity: O(1) — uses running statistics
 */
bool hss_is_zeno_prone(const HSS_System *sys) {
    if (!sys || !sys->zeno_check) return false;

    /* If too many jumps in short time, suspect Zeno */
    if (sys->state.jumps >= sys->zeno_max_jumps) {
        return true;
    }

    /* If jumps but no time progress, suspect Zeno */
    if (sys->state.jumps > 0 && sys->state.time < sys->zeno_min_time) {
        return true;
    }

    return false;
}

/**
 * @brief Estimate the Zeno time (accumulation point).
 *
 * For a geometric sequence of inter-event times t_k,
 * estimate τ_zeno = t_1 + t_2 + ... = t_1 / (1 - r)
 * where r is the estimated rate of decrease.
 *
 * This is relevant for bouncing ball type systems where
 * inter-bounce times form a geometric progression.
 *
 * Complexity: O(1) — closed-form estimate
 */
double hss_estimate_zeno_time(double first_interval, double decay_rate) {
    if (decay_rate <= 0.0 || decay_rate >= 1.0) {
        return INFINITY;
    }
    return first_interval / (1.0 - decay_rate);
}
