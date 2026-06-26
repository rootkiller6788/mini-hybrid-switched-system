#include "switched_types.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* ============================================================================
 * Vector Operations (BLAS Level 1)
 * ============================================================================ */

SwitchedVector sv_create(int n) {
    SwitchedVector v;
    v.n = n;
    v.data = (double *)calloc((size_t)n, sizeof(double));
    if (!v.data && n > 0) {
        fprintf(stderr, "sv_create: allocation failed for n=%d\n", n);
    }
    return v;
}

void sv_free(SwitchedVector *v) {
    if (!v) return;
    free(v->data);
    v->data = NULL;
    v->n = 0;
}

void sv_set(SwitchedVector *v, int i, double val) {
    if (!v || !v->data || i < 0 || i >= v->n) return;
    v->data[i] = val;
}

double sv_get(const SwitchedVector *v, int i) {
    if (!v || !v->data || i < 0 || i >= v->n) return 0.0;
    return v->data[i];
}

double sv_norm(const SwitchedVector *v) {
    if (!v || !v->data) return 0.0;
    double sum_sq = 0.0;
    for (int i = 0; i < v->n; i++) {
        sum_sq += v->data[i] * v->data[i];
    }
    return sqrt(sum_sq);
}

double sv_norm_inf(const SwitchedVector *v) {
    if (!v || !v->data) return 0.0;
    double max_val = 0.0;
    for (int i = 0; i < v->n; i++) {
        double abs_val = fabs(v->data[i]);
        if (abs_val > max_val) max_val = abs_val;
    }
    return max_val;
}

double sv_dot(const SwitchedVector *a, const SwitchedVector *b) {
    if (!a || !b || !a->data || !b->data || a->n != b->n) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < a->n; i++) {
        sum += a->data[i] * b->data[i];
    }
    return sum;
}

void sv_copy(SwitchedVector *dst, const SwitchedVector *src) {
    if (!dst || !src || !dst->data || !src->data || dst->n != src->n) return;
    memcpy(dst->data, src->data, (size_t)src->n * sizeof(double));
}

void sv_scale(SwitchedVector *v, double alpha) {
    if (!v || !v->data) return;
    for (int i = 0; i < v->n; i++) {
        v->data[i] *= alpha;
    }
}

void sv_axpy(SwitchedVector *y, double alpha, const SwitchedVector *x) {
    if (!y || !x || !y->data || !x->data || y->n != x->n) return;
    for (int i = 0; i < y->n; i++) {
        y->data[i] += alpha * x->data[i];
    }
}

void sv_print(const SwitchedVector *v, FILE *fp) {
    if (!v || !fp) return;
    fprintf(fp, "Vector[%d] = [", v->n);
    for (int i = 0; i < v->n; i++) {
        fprintf(fp, " % 8.4f", v->data ? v->data[i] : 0.0);
        if (i < v->n - 1) fprintf(fp, ",");
    }
    fprintf(fp, " ]\n");
}

/* ============================================================================
 * Matrix Operations
 * ============================================================================ */

SwitchedMatrix sm_create(int rows, int cols) {
    SwitchedMatrix M;
    M.rows = rows;
    M.cols = cols;
    M.data = (double *)calloc((size_t)(rows * cols), sizeof(double));
    if (!M.data && rows > 0 && cols > 0) {
        fprintf(stderr, "sm_create: allocation failed for %dx%d\n", rows, cols);
    }
    return M;
}

void sm_free(SwitchedMatrix *M) {
    if (!M) return;
    free(M->data);
    M->data = NULL;
    M->rows = 0;
    M->cols = 0;
}

void sm_set(SwitchedMatrix *M, int i, int j, double val) {
    if (!M || !M->data || i < 0 || i >= M->rows || j < 0 || j >= M->cols) return;
    M->data[i * M->cols + j] = val;
}

double sm_get(const SwitchedMatrix *M, int i, int j) {
    if (!M || !M->data || i < 0 || i >= M->rows || j < 0 || j >= M->cols) return 0.0;
    return M->data[i * M->cols + j];
}

void sm_identity(SwitchedMatrix *M) {
    if (!M || !M->data || M->rows != M->cols) return;
    int n = M->rows;
    memset(M->data, 0, (size_t)(n * n) * sizeof(double));
    for (int i_ = 0; i_ < n; i_++) {
        M->data[i_ * n + i_] = 1.0;
    }
}

void sm_copy(SwitchedMatrix *dst, const SwitchedMatrix *src) {
    if (!dst || !src || !dst->data || !src->data) return;
    if (dst->rows != src->rows || dst->cols != src->cols) return;
    memcpy(dst->data, src->data, (size_t)(src->rows * src->cols) * sizeof(double));
}

void sm_transpose(SwitchedMatrix *dst, const SwitchedMatrix *src) {
    if (!dst || !src || !dst->data || !src->data) return;
    if (dst->rows != src->cols || dst->cols != src->rows) return;
    for (int i = 0; i < src->rows; i++) {
        for (int j = 0; j < src->cols; j++) {
            dst->data[j * dst->cols + i] = src->data[i * src->cols + j];
        }
    }
}

void sm_add(SwitchedMatrix *C, const SwitchedMatrix *A, const SwitchedMatrix *B) {
    if (!C || !A || !B || !C->data || !A->data || !B->data) return;
    if (A->rows != B->rows || A->cols != B->cols) return;
    if (C->rows != A->rows || C->cols != A->cols) return;
    int n = A->rows * A->cols;
    for (int i = 0; i < n; i++) {
        C->data[i] = A->data[i] + B->data[i];
    }
}

void sm_sub(SwitchedMatrix *C, const SwitchedMatrix *A, const SwitchedMatrix *B) {
    if (!C || !A || !B || !C->data || !A->data || !B->data) return;
    if (A->rows != B->rows || A->cols != B->cols) return;
    if (C->rows != A->rows || C->cols != A->cols) return;
    int n = A->rows * A->cols;
    for (int i = 0; i < n; i++) {
        C->data[i] = A->data[i] - B->data[i];
    }
}

void sm_mul(SwitchedMatrix *C, const SwitchedMatrix *A, const SwitchedMatrix *B) {
    if (!C || !A || !B || !C->data || !A->data || !B->data) return;
    if (A->cols != B->rows) return;
    if (C->rows != A->rows || C->cols != B->cols) return;
    int m = A->rows, k = A->cols, n = B->cols;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int p = 0; p < k; p++) {
                sum += A->data[i * A->cols + p] * B->data[p * B->cols + j];
            }
            C->data[i * C->cols + j] = sum;
        }
    }
}

void sm_mul_scalar(SwitchedMatrix *M, double scalar) {
    if (!M || !M->data) return;
    int n = M->rows * M->cols;
    for (int i = 0; i < n; i++) {
        M->data[i] *= scalar;
    }
}

void sm_matvec_mul(SwitchedVector *y, const SwitchedMatrix *A, const SwitchedVector *x) {
    if (!y || !A || !x || !y->data || !A->data || !x->data) return;
    if (A->cols != x->n || A->rows != y->n) return;
    for (int i = 0; i < A->rows; i++) {
        double sum = 0.0;
        for (int j = 0; j < A->cols; j++) {
            sum += A->data[i * A->cols + j] * x->data[j];
        }
        y->data[i] = sum;
    }
}

double sm_trace(const SwitchedMatrix *M) {
    if (!M || !M->data || M->rows != M->cols) return 0.0;
    double tr = 0.0;
    for (int i = 0; i < M->rows; i++) {
        tr += M->data[i * M->cols + i];
    }
    return tr;
}

double sm_det_2x2(const SwitchedMatrix *M) {
    if (!M || !M->data || M->rows != 2 || M->cols != 2) return 0.0;
    return M->data[0] * M->data[3] - M->data[1] * M->data[2];
}

double sm_det_3x3(const SwitchedMatrix *M) {
    if (!M || !M->data || M->rows != 3 || M->cols != 3) return 0.0;
    double a = M->data[0], b = M->data[1], c = M->data[2];
    double d = M->data[3], e = M->data[4], f = M->data[5];
    double g = M->data[6], h = M->data[7], i = M->data[8];
    return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
}

void sm_inv_2x2(SwitchedMatrix *inv, const SwitchedMatrix *M) {
    if (!inv || !M || !inv->data || !M->data) return;
    if (inv->rows != 2 || inv->cols != 2 || M->rows != 2 || M->cols != 2) return;
    double det = sm_det_2x2(M);
    if (fabs(det) < 1e-15) return;
    double inv_det = 1.0 / det;
    inv->data[0] =  M->data[3] * inv_det;
    inv->data[1] = -M->data[1] * inv_det;
    inv->data[2] = -M->data[2] * inv_det;
    inv->data[3] =  M->data[0] * inv_det;
}

bool sm_is_symmetric(const SwitchedMatrix *M) {
    if (!M || !M->data || M->rows != M->cols) return false;
    int n = M->rows;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (fabs(M->data[i * n + j] - M->data[j * n + i]) > 1e-12) {
                return false;
            }
        }
    }
    return true;
}

bool sm_is_positive_definite(const SwitchedMatrix *M) {
    if (!M || !M->data || M->rows != M->cols) return false;
    if (!sm_is_symmetric(M)) return false;
    int n = M->rows;
    for (int k = 1; k <= n && k <= 3; k++) {
        double det_val = 0.0;
        if (k == 1) {
            det_val = M->data[0];
        } else if (k == 2) {
            det_val = M->data[0] * M->data[n + 1] - M->data[1] * M->data[n];
        } else if (k == 3 && n >= 3) {
            double a = M->data[0],          b = M->data[1],          c = M->data[2];
            double d = M->data[n],          e = M->data[n + 1],      f = M->data[n + 2];
            double g = M->data[2 * n],      h = M->data[2 * n + 1],  i = M->data[2 * n + 2];
            det_val = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
        } else {
            det_val = 1.0;
        }
        if (det_val <= 1e-12) return false;
    }
    return true;
}

double sm_frobenius_norm(const SwitchedMatrix *M) {
    if (!M || !M->data) return 0.0;
    double sum_sq = 0.0;
    int n = M->rows * M->cols;
    for (int i = 0; i < n; i++) {
        sum_sq += M->data[i] * M->data[i];
    }
    return sqrt(sum_sq);
}

void sm_print(const SwitchedMatrix *M, FILE *fp) {
    if (!M || !fp) return;
    fprintf(fp, "Matrix[%dx%d]:\n", M->rows, M->cols);
    for (int i = 0; i < M->rows; i++) {
        fprintf(fp, "  ");
        for (int j = 0; j < M->cols; j++) {
            fprintf(fp, " % 9.4f", M->data ? M->data[i * M->cols + j] : 0.0);
        }
        fprintf(fp, "\n");
    }
}

void sm_commutator(SwitchedMatrix *bracket, const SwitchedMatrix *A, const SwitchedMatrix *B) {
    if (!bracket || !A || !B) return;
    if (A->rows != A->cols || B->rows != B->cols || A->rows != B->rows) return;
    if (bracket->rows != A->rows || bracket->cols != A->cols) return;
    int n = A->rows;
    SwitchedMatrix AB = sm_create(n, n);
    SwitchedMatrix BA = sm_create(n, n);
    sm_mul(&AB, A, B);
    sm_mul(&BA, B, A);
    sm_sub(bracket, &AB, &BA);
    sm_free(&AB);
    sm_free(&BA);
}

/* ============================================================================
 * Switching Signal Operations
 * ============================================================================ */

SwitchingSignal* ssig_create(int capacity) {
    SwitchingSignal *sig = (SwitchingSignal *)malloc(sizeof(SwitchingSignal));
    if (!sig) return NULL;
    sig->type = SSIG_ARBITRARY;
    sig->capacity = capacity;
    sig->n_switches = 0;
    sig->current_mode = 0;
    sig->total_time = 0.0;
    sig->mode_sequence = (int *)malloc((size_t)capacity * sizeof(int));
    sig->switch_times = (double *)malloc((size_t)capacity * sizeof(double));
    if (!sig->mode_sequence || !sig->switch_times) {
        free(sig->mode_sequence); free(sig->switch_times); free(sig);
        return NULL;
    }
    sig->mode_sequence[0] = 0;
    sig->switch_times[0] = 0.0;
    return sig;
}

void ssig_free(SwitchingSignal *sig) {
    if (!sig) return;
    free(sig->mode_sequence);
    free(sig->switch_times);
    free(sig);
}

void ssig_record_switch(SwitchingSignal *sig, int new_mode, double time) {
    if (!sig) return;
    if (sig->n_switches + 1 >= sig->capacity) {
        int new_cap = sig->capacity * 2;
        int *nm = (int *)realloc(sig->mode_sequence, (size_t)new_cap * sizeof(int));
        double *nt = (double *)realloc(sig->switch_times, (size_t)new_cap * sizeof(double));
        if (!nm || !nt) return;
        sig->mode_sequence = nm;
        sig->switch_times = nt;
        sig->capacity = new_cap;
    }
    sig->n_switches++;
    sig->mode_sequence[sig->n_switches] = new_mode;
    sig->switch_times[sig->n_switches] = time;
    sig->current_mode = new_mode;
    sig->total_time = time;
}

void ssig_set_type(SwitchingSignal *sig, SwitchingSignalType type) {
    if (!sig) return;
    sig->type = type;
}

int ssig_active_mode_at(const SwitchingSignal *sig, double t) {
    if (!sig || sig->n_switches < 0) return -1;
    if (t <= sig->switch_times[0]) return sig->mode_sequence[0];
    int lo = 0, hi = sig->n_switches;
    while (lo < hi) {
        int mid = (lo + hi + 1) / 2;
        if (sig->switch_times[mid] <= t) lo = mid;
        else hi = mid - 1;
    }
    return sig->mode_sequence[lo];
}

void ssig_reset(SwitchingSignal *sig) {
    if (!sig) return;
    sig->n_switches = 0;
    sig->current_mode = 0;
    sig->total_time = 0.0;
    if (sig->mode_sequence) sig->mode_sequence[0] = 0;
    if (sig->switch_times) sig->switch_times[0] = 0.0;
}

void ssig_print(const SwitchingSignal *sig, FILE *fp) {
    if (!sig || !fp) return;
    fprintf(fp, "SwitchingSignal: type=%d, switches=%d\n", sig->type, sig->n_switches);
    for (int k = 0; k <= sig->n_switches && k < sig->capacity; k++)
        fprintf(fp, "  t=%.4f: mode=%d\n", sig->switch_times[k], sig->mode_sequence[k]);
}

/* ============================================================================
 * Subsystem Operations
 * ============================================================================ */

SwitchedSubsystem* ssub_create(int mode_id, int state_dim) {
    SwitchedSubsystem *sub = (SwitchedSubsystem *)malloc(sizeof(SwitchedSubsystem));
    if (!sub) return NULL;
    sub->mode_id = mode_id;
    sub->A = sm_create(state_dim, state_dim);
    sub->B = sm_create(state_dim, 1);
    sub->spectral_radius = 0.0;
    sub->is_hurwitz = false;
    sub->is_schur = false;
    snprintf(sub->description, sizeof(sub->description), "Subsystem %d", mode_id);
    return sub;
}

void ssub_free(SwitchedSubsystem *sub) {
    if (!sub) return;
    sm_free(&sub->A);
    sm_free(&sub->B);
    free(sub);
}

void ssub_set_hurwitz_check(SwitchedSubsystem *sub) {
    if (!sub) return;
    sub->spectral_radius = spectral_radius(&sub->A);
    sub->is_hurwitz = is_hurwitz_matrix(&sub->A);
    sub->is_schur = (sub->spectral_radius < 1.0);
}

void ssub_print(const SwitchedSubsystem *sub, FILE *fp) {
    if (!sub || !fp) return;
    fprintf(fp, "Subsystem %d: %s\n", sub->mode_id, sub->description);
    fprintf(fp, "  Hurwitz=%s, rho=%.6f\n", sub->is_hurwitz ? "yes" : "no", sub->spectral_radius);
    sm_print(&sub->A, fp);
}

/* ============================================================================
 * Switched System Operations
 * ============================================================================ */

SwitchedSystem* sss_create(const char *name, int state_dim, int n_modes) {
    SwitchedSystem *sys = (SwitchedSystem *)malloc(sizeof(SwitchedSystem));
    if (!sys) return NULL;
    strncpy(sys->name, name, sizeof(sys->name) - 1);
    sys->name[sizeof(sys->name) - 1] = '\0';
    sys->state_dim = state_dim;
    sys->n_modes = n_modes;
    sys->t = 0.0;
    sys->modes = (SwitchedSubsystem **)calloc((size_t)n_modes, sizeof(SwitchedSubsystem *));
    if (!sys->modes) { free(sys); return NULL; }
    sys->signal = NULL;
    sys->state = sv_create(state_dim);
    sys->stability = SSTAB_UNSTABLE;
    sys->eigenvalue_dim = state_dim;
    sys->eigenvalues = (double *)calloc((size_t)(state_dim * 2), sizeof(double));
    sys->lf_type = CLF_COMMON;
    sys->lyap_matrices = NULL;
    sys->lyap_values = (double *)calloc((size_t)n_modes, sizeof(double));
    sys->mode_active = (bool *)calloc((size_t)n_modes, sizeof(bool));
    sys->avg_dwell_time = 0.0;
    sys->min_dwell_time = 0.0;
    sys->dwell_samples = NULL;
    sys->dwell_count = 0;
    return sys;
}

void sss_free(SwitchedSystem *sys) {
    if (!sys) return;
    if (sys->modes) {
        for (int i = 0; i < sys->n_modes; i++) ssub_free(sys->modes[i]);
        free(sys->modes);
    }
    sv_free(&sys->state);
    free(sys->eigenvalues);
    free(sys->lyap_values);
    free(sys->mode_active);
    free(sys->dwell_samples);
    if (sys->lyap_matrices) {
        for (int i = 0; i < sys->n_modes; i++) {
            if (sys->lyap_matrices[i]) { sm_free(sys->lyap_matrices[i]); free(sys->lyap_matrices[i]); }
        }
        free(sys->lyap_matrices);
    }
    free(sys);
}

void sss_add_subsystem(SwitchedSystem *sys, int mode_id, const SwitchedMatrix *A) {
    if (!sys || !A || mode_id < 0 || mode_id >= sys->n_modes) return;
    SwitchedSubsystem *sub = ssub_create(mode_id, sys->state_dim);
    if (!sub) return;
    sm_copy(&sub->A, A);
    ssub_set_hurwitz_check(sub);
    sys->modes[mode_id] = sub;
}

void sss_set_signal(SwitchedSystem *sys, SwitchingSignal *sig) {
    if (sys) sys->signal = sig;
}

void sss_set_initial_state(SwitchedSystem *sys, const SwitchedVector *x0) {
    if (!sys || !x0 || x0->n != sys->state_dim) return;
    sv_copy(&sys->state, x0);
    sys->t = 0.0;
}

void sss_step(SwitchedSystem *sys, double dt) {
    if (!sys || !sys->signal || !sys->modes) return;
    int mode = ssig_active_mode_at(sys->signal, sys->t);
    if (mode < 0 || mode >= sys->n_modes) return;
    SwitchedVector Ax = sv_create(sys->state_dim);
    sm_matvec_mul(&Ax, &sys->modes[mode]->A, &sys->state);
    sv_axpy(&sys->state, dt, &Ax);
    sv_free(&Ax);
    sys->t += dt;
}

void sss_simulate(SwitchedSystem *sys, const SolverConfig *cfg, const SwitchSequence *seq) {
    if (!sys || !cfg || !seq) return;
    if (!sys->signal) sys->signal = ssig_create(seq->length + 2);
    double t = 0.0;
    int si = 0;
    double tim = 0.0;
    sys->signal->mode_sequence[0] = seq->mode_order[0];
    sys->signal->switch_times[0] = 0.0;
    sys->signal->n_switches = 0;
    sys->signal->current_mode = seq->mode_order[0];
    for (int step = 0; step < cfg->max_steps && t < cfg->t_end; step++) {
        double dt = cfg->dt;
        if (t + dt > cfg->t_end) dt = cfg->t_end - t;
        double rem = seq->durations[si] - tim;
        if (dt > rem && rem > cfg->event_tol) {
            sss_step(sys, rem);
            t += rem;
            si++;
            if (si >= seq->length) { if (seq->repeat) si = 0; else break; }
            ssig_record_switch(sys->signal, seq->mode_order[si], t);
            tim = 0.0;
            double dtr = dt - rem;
            sss_step(sys, dtr);
            t += dtr;
            tim += dtr;
        } else {
            sss_step(sys, dt);
            t += dt;
            tim += dt;
            if (fabs(tim - seq->durations[si]) < cfg->event_tol) {
                si++;
                if (si < seq->length) {
                    ssig_record_switch(sys->signal, seq->mode_order[si], t);
                    tim = 0.0;
                } else if (seq->repeat) {
                    si = 0;
                    ssig_record_switch(sys->signal, seq->mode_order[si], t);
                    tim = 0.0;
                }
            }
        }
    }
    sys->signal->total_time = t;
}

void sss_update_lyapunov_values(SwitchedSystem *sys) {
    if (!sys || !sys->lyap_matrices || !sys->lyap_values) return;
    for (int i = 0; i < sys->n_modes; i++) {
        if (sys->lyap_matrices[i]) {
            SwitchedVector Px = sv_create(sys->state_dim);
            sm_matvec_mul(&Px, sys->lyap_matrices[i], &sys->state);
            sys->lyap_values[i] = sv_dot(&sys->state, &Px);
            sv_free(&Px);
        }
    }
}

void sss_print(const SwitchedSystem *sys, FILE *fp) {
    if (!sys || !fp) return;
    fprintf(fp, "SwitchedSystem: %s\n", sys->name);
    fprintf(fp, "  dim=%d, modes=%d, t=%.6f, stability=%d\n",
            sys->state_dim, sys->n_modes, sys->t, sys->stability);
    sv_print(&sys->state, fp);
    for (int i = 0; i < sys->n_modes; i++) ssub_print(sys->modes[i], fp);
}

/* ============================================================================
 * Lie Algebra Operations
 * ============================================================================ */

LieAlgebraCondition* lie_check_create(int n_modes, int dim) {
    LieAlgebraCondition *la = (LieAlgebraCondition *)malloc(sizeof(LieAlgebraCondition));
    if (!la) return NULL;
    la->is_solvable = false;
    la->is_nilpotent = false;
    la->pair_commute = true;
    la->simultaneously_triang = false;
    la->n_modes = n_modes;
    la->dim = dim;
    la->lie_brackets = (double *)calloc((size_t)(n_modes * n_modes * dim * dim), sizeof(double));
    return la;
}

void lie_check_free(LieAlgebraCondition *la) {
    if (!la) return;
    free(la->lie_brackets);
    free(la);
}

void lie_compute_bracket(SwitchedMatrix *bracket, const SwitchedMatrix *A, const SwitchedMatrix *B) {
    sm_commutator(bracket, A, B);
}

bool lie_check_commute(const SwitchedMatrix *A, const SwitchedMatrix *B) {
    if (!A || !B || A->rows != B->rows) return false;
    int n = A->rows;
    SwitchedMatrix bracket = sm_create(n, n);
    sm_commutator(&bracket, A, B);
    double fnorm = sm_frobenius_norm(&bracket);
    sm_free(&bracket);
    return (fnorm < 1e-10);
}

bool lie_check_simult_triang(SwitchedMatrix **A, int n_modes, int dim) {
    (void)dim; /* reserved for future dimension-dependent checks */
    if (n_modes <= 1) return true;
    for (int i = 0; i < n_modes; i++)
        for (int j = i + 1; j < n_modes; j++)
            if (!lie_check_commute(A[i], A[j])) return false;
    return true;
}

bool lie_condition_solvable(LieAlgebraCondition *la, SwitchedMatrix **A, int n_modes, int dim) {
    if (!la || !A) return false;
    la->pair_commute = true;
    for (int i = 0; i < n_modes && la->pair_commute; i++)
        for (int j = i + 1; j < n_modes && la->pair_commute; j++)
            if (!lie_check_commute(A[i], A[j])) la->pair_commute = false;
    if (la->pair_commute) {
        la->is_solvable = true;
        la->is_nilpotent = true;
        la->simultaneously_triang = true;
        return true;
    }
    la->simultaneously_triang = lie_check_simult_triang(A, n_modes, dim);
    la->is_solvable = la->simultaneously_triang;
    return la->is_solvable;
}

/* ============================================================================
 * QR Algorithm for Eigenvalues
 * ============================================================================ */

QRWorkspace* qr_create(int n) {
    QRWorkspace *qr = (QRWorkspace *)malloc(sizeof(QRWorkspace));
    if (!qr) return NULL;
    qr->H = sm_create(n, n);
    qr->Q = sm_create(n, n);
    qr->wr = (double *)calloc((size_t)n, sizeof(double));
    qr->wi = (double *)calloc((size_t)n, sizeof(double));
    qr->n = n;
    qr->iters = 0;
    qr->converged = false;
    return qr;
}

void qr_free(QRWorkspace *qr) {
    if (!qr) return;
    sm_free(&qr->H);
    sm_free(&qr->Q);
    free(qr->wr);
    free(qr->wi);
    free(qr);
}

int qr_eigenvalues(QRWorkspace *qr, SwitchedMatrix *M, EigenvalueResult *results) {
    if (!qr || !M || !results) return -1;
    if (M->rows != M->cols || M->rows != qr->n) return -1;
    int n = qr->n;
    sm_copy(&qr->H, M);
    qr->iters = 0;
    qr->converged = false;

    /* n == 1: trivial eigenvalue */
    if (n == 1) {
        results[0].real = qr->H.data[0];
        results[0].imag = 0.0;
        results[0].magnitude = fabs(results[0].real);
        results[0].stable = (results[0].real < 0.0);
        qr->converged = true;
        return 1;
    }

    /* n == 2: quadratic formula for eigenvalues */
    if (n == 2) {
        double a = qr->H.data[0], b = qr->H.data[1];
        double c = qr->H.data[2], d = qr->H.data[3];
        double tr = a + d;
        double det = a * d - b * c;
        double disc = tr * tr - 4.0 * det;
        if (disc >= 0.0) {
            double sqrt_d = sqrt(disc);
            results[0].real = (tr + sqrt_d) / 2.0; results[0].imag = 0.0;
            results[1].real = (tr - sqrt_d) / 2.0; results[1].imag = 0.0;
        } else {
            double re = tr / 2.0, im = sqrt(-disc) / 2.0;
            results[0].real = re; results[0].imag = im;
            results[1].real = re; results[1].imag = -im;
        }
        for (int i = 0; i < 2; i++) {
            results[i].magnitude = sqrt(results[i].real * results[i].real +
                                        results[i].imag * results[i].imag);
            results[i].stable = (results[i].real < 0.0);
        }
        qr->converged = true;
        return 2;
    }

    /* n >= 3: QR iteration with Wilkinson shift */
    int max_iter = QR_MAX_ITERS;
    int m = n;

    while (m > 2 && qr->iters < max_iter) {
        double hm = qr->H.data[(m-1) * n + (m-1)];
        double h1 = qr->H.data[(m-2) * n + (m-2)];
        double hm1 = qr->H.data[(m-1) * n + (m-2)];
        double h1m = qr->H.data[(m-2) * n + (m-1)];
        double tr2 = h1 + hm;
        double det2 = h1 * hm - h1m * hm1;
        double disc2 = tr2 * tr2 - 4.0 * det2;
        double mu;

        if (disc2 >= 0.0) {
            double s = sqrt(disc2);
            double mu1 = (tr2 + s) / 2.0;
            double mu2 = (tr2 - s) / 2.0;
            mu = (fabs(mu1 - hm) < fabs(mu2 - hm)) ? mu1 : mu2;
        } else {
            mu = tr2 / 2.0;
        }

        /* Apply shift */
        for (int i = 0; i < m; i++) qr->H.data[i * n + i] -= mu;

        /* QR via Householder */
        for (int k = 0; k < m - 1; k++) {
            double xn = 0.0;
            for (int i = k; i < m; i++) {
                double val = qr->H.data[i * n + k];
                xn += val * val;
            }
            xn = sqrt(xn);
            if (xn < QR_TOL) continue;

            double alpha = qr->H.data[k * n + k];
            if (alpha > 0) xn = -xn;

            double *v = (double *)malloc((size_t)m * sizeof(double));
            v[k] = qr->H.data[k * n + k] - xn;
            double vns = v[k] * v[k];
            for (int i = k + 1; i < m; i++) {
                v[i] = qr->H.data[i * n + k];
                vns += v[i] * v[i];
            }
            double beta = 2.0 / vns;

            for (int j = k; j < m; j++) {
                double vdc = 0.0;
                for (int i = k; i < m; i++) vdc += v[i] * qr->H.data[i * n + j];
                double fac = beta * vdc;
                for (int i = k; i < m; i++) qr->H.data[i * n + j] -= fac * v[i];
            }
            for (int ii = 0; ii < m; ii++) {
                double vdr = 0.0;
                for (int jj = k; jj < m; jj++) vdr += v[jj] * qr->H.data[ii * n + jj];
                double fac = beta * vdr;
                for (int jj = k; jj < m; jj++) qr->H.data[ii * n + jj] -= fac * v[jj];
            }
            free(v);
        }

        /* Restore shift */
        for (int i = 0; i < m; i++) qr->H.data[i * n + i] += mu;
        qr->iters++;

        /* Deflation check */
        if (fabs(qr->H.data[(m-1) * n + (m-2)]) < QR_TOL) {
            int idx = m - 1;
            results[idx].real = qr->H.data[idx * n + idx];
            results[idx].imag = 0.0;
            results[idx].magnitude = fabs(results[idx].real);
            results[idx].stable = (results[idx].real < 0.0);
            m--;
        }
    }

    /* Handle remaining 2x2 or 1x1 block */
    if (m == 2) {
        double a = qr->H.data[0], b = qr->H.data[1];
        double c = qr->H.data[n], d = qr->H.data[n + 1];
        double tr = a + d, det = a * d - b * c, disc = tr * tr - 4.0 * det;
        if (disc >= 0.0) {
            double s = sqrt(disc);
            results[0].real = (tr + s) / 2.0; results[0].imag = 0.0;
            results[1].real = (tr - s) / 2.0; results[1].imag = 0.0;
        } else {
            double re = tr / 2.0, im = sqrt(-disc) / 2.0;
            results[0].real = re; results[0].imag = im;
            results[1].real = re; results[1].imag = -im;
        }
        for (int i = 0; i < 2; i++) {
            results[i].magnitude = sqrt(results[i].real * results[i].real +
                                        results[i].imag * results[i].imag);
            results[i].stable = (results[i].real < 0.0);
        }
    } else if (m == 1) {
        results[0].real = qr->H.data[0];
        results[0].imag = 0.0;
        results[0].magnitude = fabs(results[0].real);
        results[0].stable = (results[0].real < 0.0);
    }

    qr->converged = (qr->iters < max_iter);
    for (int i = 0; i < n; i++) {
        qr->wr[i] = results[i].real;
        qr->wi[i] = results[i].imag;
    }
    return n;
}

double spectral_radius(const SwitchedMatrix *M) {
    if (!M || M->rows != M->cols) return 0.0;
    int n = M->rows;
    QRWorkspace *qr = qr_create(n);
    EigenvalueResult *eig = (EigenvalueResult *)malloc((size_t)n * sizeof(EigenvalueResult));
    qr_eigenvalues(qr, (SwitchedMatrix *)M, eig);
    double rho = 0.0;
    for (int i = 0; i < n; i++) {
        if (eig[i].magnitude > rho) rho = eig[i].magnitude;
    }
    free(eig);
    qr_free(qr);
    return rho;
}

bool is_hurwitz_matrix(const SwitchedMatrix *M) {
    if (!M || M->rows != M->cols) return false;
    int n = M->rows;
    QRWorkspace *qr = qr_create(n);
    EigenvalueResult *eig = (EigenvalueResult *)malloc((size_t)n * sizeof(EigenvalueResult));
    qr_eigenvalues(qr, (SwitchedMatrix *)M, eig);
    bool hurwitz = true;
    for (int i = 0; i < n; i++) {
        if (eig[i].real >= -QR_TOL) { hurwitz = false; break; }
    }
    free(eig);
    qr_free(qr);
    return hurwitz;
}
