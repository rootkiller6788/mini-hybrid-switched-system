#include "etc_core.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ETC_HISTORY_DEFAULT_CAP 1024

/* ============================================================================
 * Internal Helpers
 * ============================================================================ */

static double etc_max(double a, double b) { return a > b ? a : b; }

static double etc_abs(double x) { return x < 0.0 ? -x : x; }

/* Internal: default quadratic trigger function */
static double etc_default_trigger(const ETCVector* x, const ETCVector* e,
                                   double sigma, double epsilon) {
    (void)epsilon;
    double x_norm_sq = etc_vector_dot(x, x);
    double e_norm_sq = etc_vector_dot(e, e);
    return e_norm_sq - sigma * x_norm_sq;
}

/* Internal: compute maximum absolute eigenvalue for 2x2 matrix */
static void etc_eigen_2x2_direct(double a, double b, double c, double d,
                                  double* re, double* im) {
    double trace = a + d;
    double det = a * d - b * c;
    double disc = trace * trace - 4.0 * det;
    if (disc >= 0.0) {
        double sqrt_disc = sqrt(disc);
        re[0] = 0.5 * (trace + sqrt_disc);
        im[0] = 0.0;
        re[1] = 0.5 * (trace - sqrt_disc);
        im[1] = 0.0;
    } else {
        re[0] = 0.5 * trace;
        im[0] = 0.5 * sqrt(-disc);
        re[1] = re[0];
        im[1] = -im[0];
    }
}

/* Internal: power iteration for dominant eigenvalue */
static double etc_power_iteration(const ETCMatrix* A, int max_iter, double tol) {
    int n = A->rows;
    double* v = (double*)malloc(n * sizeof(double));
    double* Av = (double*)malloc(n * sizeof(double));
    if (!v || !Av) { free(v); free(Av); return 0.0; }
    for (int i = 0; i < n; i++) v[i] = 1.0 / sqrt((double)n);
    double lambda = 0.0, lambda_old = 0.0;
    for (int iter = 0; iter < max_iter; iter++) {
        double v_norm = 0.0;
        for (int i = 0; i < n; i++) {
            Av[i] = 0.0;
            for (int j = 0; j < n; j++)
                Av[i] += A->data[i * n + j] * v[j];
        }
        lambda = 0.0;
        for (int i = 0; i < n; i++) {
            lambda += v[i] * Av[i];
            v_norm += Av[i] * Av[i];
        }
        v_norm = sqrt(v_norm);
        if (v_norm < 1e-15) break;
        double inv_norm = 1.0 / v_norm;
        for (int i = 0; i < n; i++) v[i] = Av[i] * inv_norm;
        if (etc_abs(lambda - lambda_old) < tol * etc_max(1.0, etc_abs(lambda)))
            break;
        lambda_old = lambda;
    }
    free(v); free(Av);
    return lambda;
}

/* ============================================================================
 * Vector Operations
 * ============================================================================ */

ETCVector etc_vector_create(int dim) {
    ETCVector v;
    v.dim = dim;
    v.data = (double*)calloc((size_t)dim, sizeof(double));
    return v;
}

void etc_vector_free(ETCVector* v) {
    if (v && v->data) { free(v->data); v->data = NULL; v->dim = 0; }
}

double etc_vector_norm(const ETCVector* v) {
    if (!v || !v->data || v->dim <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < v->dim; i++) sum += v->data[i] * v->data[i];
    return sqrt(sum);
}

double etc_vector_dot(const ETCVector* v1, const ETCVector* v2) {
    if (!v1 || !v2 || !v1->data || !v2->data || v1->dim != v2->dim) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < v1->dim; i++) sum += v1->data[i] * v2->data[i];
    return sum;
}

void etc_vector_sub(const ETCVector* v1, const ETCVector* v2, ETCVector* v_out) {
    if (!v1 || !v2 || !v_out || v1->dim != v2->dim || v1->dim != v_out->dim
        || !v1->data || !v2->data || !v_out->data) return;
    for (int i = 0; i < v1->dim; i++)
        v_out->data[i] = v1->data[i] - v2->data[i];
}

void etc_vector_scale(double s, const ETCVector* v, ETCVector* v_out) {
    if (!v || !v_out || !v->data || !v_out->data || v->dim != v_out->dim) return;
    for (int i = 0; i < v->dim; i++)
        v_out->data[i] = s * v->data[i];
}

void etc_vector_add(const ETCVector* v1, const ETCVector* v2, ETCVector* v_out) {
    if (!v1 || !v2 || !v_out || v1->dim != v2->dim || v1->dim != v_out->dim
        || !v1->data || !v2->data || !v_out->data) return;
    for (int i = 0; i < v1->dim; i++)
        v_out->data[i] = v1->data[i] + v2->data[i];
}

void etc_vector_copy(const ETCVector* src, ETCVector* dst) {
    if (!src || !dst || !src->data || !dst->data || src->dim != dst->dim) return;
    for (int i = 0; i < src->dim; i++) dst->data[i] = src->data[i];
}

/* ============================================================================
 * Matrix Operations
 * ============================================================================ */

ETCMatrix etc_matrix_create(int rows, int cols) {
    ETCMatrix m;
    m.rows = rows;
    m.cols = cols;
    m.data = (double*)calloc((size_t)(rows * cols), sizeof(double));
    return m;
}

void etc_matrix_free(ETCMatrix* m) {
    if (m && m->data) { free(m->data); m->data = NULL; m->rows = 0; m->cols = 0; }
}

void etc_matrix_vec_mul(const ETCMatrix* A, const ETCVector* x, ETCVector* y) {
    if (!A || !x || !y || !A->data || !x->data || !y->data) return;
    if (A->cols != x->dim || A->rows != y->dim) return;
    for (int i = 0; i < A->rows; i++) {
        y->data[i] = 0.0;
        for (int j = 0; j < A->cols; j++)
            y->data[i] += A->data[i * A->cols + j] * x->data[j];
    }
}

void etc_matrix_mul(const ETCMatrix* A, const ETCMatrix* B, ETCMatrix* C) {
    if (!A || !B || !C || !A->data || !B->data || !C->data) return;
    if (A->cols != B->rows || C->rows != A->rows || C->cols != B->cols) return;
    int m = A->rows, k = A->cols, n = B->cols;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int p = 0; p < k; p++)
                sum += A->data[i * k + p] * B->data[p * n + j];
            C->data[i * n + j] = sum;
        }
    }
}

void etc_matrix_add(const ETCMatrix* A, const ETCMatrix* B, ETCMatrix* C) {
    if (!A || !B || !C || !A->data || !B->data || !C->data) return;
    if (A->rows != B->rows || A->cols != B->cols
        || C->rows != A->rows || C->cols != A->cols) return;
    for (int i = 0; i < A->rows * A->cols; i++)
        C->data[i] = A->data[i] + B->data[i];
}

void etc_matrix_scale(double s, const ETCMatrix* A, ETCMatrix* B) {
    if (!A || !B || !A->data || !B->data) return;
    if (A->rows != B->rows || A->cols != B->cols) return;
    for (int i = 0; i < A->rows * A->cols; i++)
        B->data[i] = s * A->data[i];
}

void etc_matrix_eigenvalues_2x2(const ETCMatrix* A, double* re, double* im) {
    if (!A || !re || !im || A->rows != 2 || A->cols != 2) return;
    etc_eigen_2x2_direct(A->data[0], A->data[1], A->data[2], A->data[3], re, im);
}

double etc_matrix_spectral_radius(const ETCMatrix* A) {
    if (!A || !A->data || A->rows != A->cols || A->rows <= 0) return 0.0;
    int n = A->rows;
    if (n == 1) return etc_abs(A->data[0]);
    if (n == 2) {
        double re[2], im[2];
        etc_matrix_eigenvalues_2x2(A, re, im);
        double r1 = sqrt(re[0]*re[0] + im[0]*im[0]);
        double r2 = sqrt(re[1]*re[1] + im[1]*im[1]);
        return etc_max(r1, r2);
    }
    return etc_abs(etc_power_iteration(A, 200, 1e-10));
}

bool etc_matrix_is_symmetric(const ETCMatrix* A) {
    if (!A || A->rows != A->cols) return false;
    int n = A->rows;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (etc_abs(A->data[i * n + j] - A->data[j * n + i]) > 1e-12)
                return false;
    return true;
}

bool etc_matrix_is_positive_definite(const ETCMatrix* A) {
    if (!A || !A->data || A->rows != A->cols || A->rows <= 0) return false;
    if (!etc_matrix_is_symmetric(A)) return false;
    int n = A->rows;
    /* Check leading principal minors for small matrices */
    if (n <= 8) {
        for (int k = 1; k <= n; k++) {
            /* Compute k×k leading principal minor via simple elimination */
            double* L = (double*)malloc((size_t)(k * k) * sizeof(double));
            double det = 1.0;
            bool ok = true;
            for (int i = 0; i < k; i++)
                for (int j = 0; j < k; j++)
                    L[i * k + j] = A->data[i * n + j];
            for (int i = 0; i < k && ok; i++) {
                if (L[i * k + i] <= 0.0) { ok = false; break; }
                det *= L[i * k + i];
                for (int j = i + 1; j < k; j++) {
                    double factor = L[j * k + i] / L[i * k + i];
                    for (int c = i; c < k; c++)
                        L[j * k + c] -= factor * L[i * k + c];
                }
            }
            free(L);
            if (!ok || det <= 0.0) return false;
        }
        return true;
    }
    /* For larger matrices: use Gershgorin circle theorem as necessary condition */
    for (int i = 0; i < n; i++) {
        double radius = 0.0;
        for (int j = 0; j < n; j++)
            if (i != j) radius += etc_abs(A->data[i * n + j]);
        if (A->data[i * n + i] <= radius) return false; /* Not strictly diagonally dominant */
    }
    /* Passed Gershgorin; also check all diagonal entries positive */
    for (int i = 0; i < n; i++)
        if (A->data[i * n + i] <= 0.0) return false;
    return true;
}

double etc_matrix_norm_2(const ETCMatrix* A) {
    if (!A || !A->data) return 0.0;
    int n = A->rows, m = A->cols;
    /* Compute A^T A spectral radius */
    ETCMatrix AtA = etc_matrix_create(m, m);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++)
                sum += A->data[k * m + i] * A->data[k * m + j];
            AtA.data[i * m + j] = sum;
        }
    }
    double rho = etc_matrix_spectral_radius(&AtA);
    etc_matrix_free(&AtA);
    return sqrt(rho);
}

/* ============================================================================
 * ETC System Core
 * ============================================================================ */

ETCSystem* etc_system_create(const double* A, const double* B,
                              const double* K, int n, int m) {
    if (!A || !B || !K || n <= 0 || m <= 0) return NULL;
    ETCSystem* sys = (ETCSystem*)calloc(1, sizeof(ETCSystem));
    if (!sys) return NULL;

    sys->n_states = n;
    sys->n_inputs = m;
    sys->A = etc_matrix_create(n, n);
    sys->B = etc_matrix_create(n, m);
    sys->K = etc_matrix_create(m, n);
    sys->Acl = etc_matrix_create(n, n);

    /* Copy matrices */
    for (int i = 0; i < n * n; i++) sys->A.data[i] = A[i];
    for (int i = 0; i < n * m; i++) sys->B.data[i] = B[i];
    for (int i = 0; i < m * n; i++) sys->K.data[i] = K[i];

    /* Precompute Acl = A + B*K */
    ETCMatrix BK = etc_matrix_create(n, n);
    etc_matrix_mul(&sys->B, &sys->K, &BK);
    etc_matrix_add(&sys->A, &BK, &sys->Acl);
    etc_matrix_free(&BK);

    /* Allocate state vectors */
    sys->x = etc_vector_create(n);
    sys->x_hat = etc_vector_create(n);
    sys->e = etc_vector_create(n);
    sys->u = etc_vector_create(m);

    /* Default trigger: quadratic */
    sys->trigger_type = ETC_TRIGGER_STATIC;
    sys->sigma = 0.1;
    sys->epsilon = 0.0;
    sys->trigger_fn = etc_default_trigger;

    /* Dynamic trigger defaults */
    sys->eta = 0.0;
    sys->eta_dot = 0.0;
    sys->theta = 0.5;

    /* Default Lyapunov */
    sys->V.P = etc_matrix_create(n, n);
    sys->V.P_norm = 0.0;
    sys->V.lambda_min_P = 0.0;
    sys->V.lambda_max_P = 0.0;
    sys->V.is_positive_definite = false;

    /* Timing */
    sys->t = 0.0;
    sys->t_last_event = 0.0;
    sys->dt = 0.001;
    sys->t_max = 10.0;
    sys->event_count = 0;
    sys->min_iet_observed = 1e100;
    sys->regime = ETC_REGIME_UNKNOWN;
    sys->zeno_bound = INFINITY;

    /* History */
    etc_history_init(&sys->history, ETC_HISTORY_DEFAULT_CAP);

    return sys;
}

void etc_system_free(ETCSystem* sys) {
    if (!sys) return;
    etc_matrix_free(&sys->A);
    etc_matrix_free(&sys->B);
    etc_matrix_free(&sys->K);
    etc_matrix_free(&sys->Acl);
    etc_matrix_free(&sys->V.P);
    etc_vector_free(&sys->x);
    etc_vector_free(&sys->x_hat);
    etc_vector_free(&sys->e);
    etc_vector_free(&sys->u);
    etc_history_free(&sys->history);
    free(sys);
}

void etc_system_set_initial_state(ETCSystem* sys, const double* x0) {
    if (!sys || !x0) return;
    for (int i = 0; i < sys->n_states; i++) {
        sys->x.data[i] = x0[i];
        sys->x_hat.data[i] = x0[i];
        sys->e.data[i] = 0.0;
    }
    sys->t = 0.0;
    sys->t_last_event = 0.0;
    sys->event_count = 0;
    sys->min_iet_observed = 1e100;
}

void etc_system_set_trigger(ETCSystem* sys, ETCTriggerType type,
                             double sigma, double epsilon,
                             double (*trigger_fn)(const ETCVector*, const ETCVector*,
                                                   double, double)) {
    if (!sys) return;
    sys->trigger_type = type;
    sys->sigma = sigma;
    sys->epsilon = epsilon;
    if (trigger_fn) sys->trigger_fn = trigger_fn;
}

void etc_system_set_dynamic_trigger(ETCSystem* sys, double eta0, double theta) {
    if (!sys) return;
    sys->eta = eta0 >= 0.0 ? eta0 : 0.0;
    sys->theta = theta > 0.0 ? theta : 0.5;
    sys->eta_dot = 0.0;
}

double etc_system_error_norm(const ETCSystem* sys) {
    if (!sys) return 0.0;
    return etc_vector_norm(&sys->e);
}

double etc_system_state_norm(const ETCSystem* sys) {
    if (!sys) return 0.0;
    return etc_vector_norm(&sys->x);
}

double etc_system_eval_trigger(const ETCSystem* sys) {
    if (!sys || !sys->trigger_fn) return -1.0;
    if (sys->trigger_type == ETC_TRIGGER_DYNAMIC) {
        double static_part = etc_default_trigger(&sys->x, &sys->e, sys->sigma, 0.0);
        return sys->theta * sys->eta + static_part;
    }
    return sys->trigger_fn(&sys->x, &sys->e, sys->sigma, sys->epsilon);
}

bool etc_system_should_trigger(const ETCSystem* sys) {
    return etc_system_eval_trigger(sys) >= 0.0;
}

void etc_system_trigger_event(ETCSystem* sys) {
    if (!sys) return;
    /* Record inter-event time */
    double iet = sys->t - sys->t_last_event;
    if (iet < sys->min_iet_observed) sys->min_iet_observed = iet;

    /* Build event record */
    ETCEvent ev;
    ev.time = sys->t;
    ev.event_index = sys->event_count;
    ev.state = etc_vector_create(sys->n_states);
    ev.control = etc_vector_create(sys->n_inputs);
    ev.error_at_event = etc_vector_create(sys->n_states);
    etc_vector_copy(&sys->x, &ev.state);
    etc_vector_copy(&sys->u, &ev.control);
    etc_vector_copy(&sys->e, &ev.error_at_event);
    ev.trigger_value = etc_system_eval_trigger(sys);
    ev.inter_event_time = iet;

    /* Append to history */
    etc_history_append(&sys->history, &ev);

    /* Reset sampled state and error */
    etc_vector_copy(&sys->x, &sys->x_hat);
    for (int i = 0; i < sys->n_states; i++) sys->e.data[i] = 0.0;

    sys->t_last_event = sys->t;
    sys->event_count++;

    /* Update control */
    etc_system_compute_control(sys);

    /* Clean event vectors */
    etc_vector_free(&ev.state);
    etc_vector_free(&ev.control);
    etc_vector_free(&ev.error_at_event);
}

void etc_system_compute_control(ETCSystem* sys) {
    if (!sys) return;
    etc_matrix_vec_mul(&sys->K, &sys->x_hat, &sys->u);
}

void etc_system_compute_derivative(const ETCSystem* sys, double* dx) {
    if (!sys || !dx) return;
    /* dx = A*x + B*u */
    ETCVector Ax = etc_vector_create(sys->n_states);
    ETCVector Bu = etc_vector_create(sys->n_states);
    etc_matrix_vec_mul(&sys->A, &sys->x, &Ax);
    etc_matrix_vec_mul(&sys->B, &sys->u, &Bu);
    for (int i = 0; i < sys->n_states; i++) dx[i] = Ax.data[i] + Bu.data[i];
    etc_vector_free(&Ax);
    etc_vector_free(&Bu);
}

bool etc_system_step_rk4(ETCSystem* sys) {
    if (!sys) return false;
    int n = sys->n_states;
    double dt = sys->dt;

    /* Check Zeno */
    if (sys->min_iet_observed < 1e-12 && sys->event_count > 100) {
        sys->regime = ETC_REGIME_ZENO;
        sys->zeno_bound = sys->t;
        return false;
    }

    /* RK4 integration */
    double* k1 = (double*)malloc((size_t)(4 * n) * sizeof(double));
    double* k2 = k1 + n;
    double* k3 = k2 + n;
    double* k4 = k3 + n;
    double* x_temp = (double*)malloc((size_t)n * sizeof(double));
    if (!k1 || !x_temp) { free(k1); free(x_temp); return false; }

    /* Save current state */
    double* x_save = (double*)malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) x_save[i] = sys->x.data[i];

    /* k1 */
    etc_system_compute_derivative(sys, k1);

    /* k2 */
    for (int i = 0; i < n; i++) sys->x.data[i] = x_save[i] + 0.5 * dt * k1[i];
    etc_system_compute_derivative(sys, k2);

    /* k3 */
    for (int i = 0; i < n; i++) sys->x.data[i] = x_save[i] + 0.5 * dt * k2[i];
    etc_system_compute_derivative(sys, k3);

    /* k4 */
    for (int i = 0; i < n; i++) sys->x.data[i] = x_save[i] + dt * k3[i];
    etc_system_compute_derivative(sys, k4);

    /* Update state */
    for (int i = 0; i < n; i++)
        sys->x.data[i] = x_save[i] + (dt / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);

    /* Update error */
    for (int i = 0; i < n; i++)
        sys->e.data[i] = sys->x_hat.data[i] - sys->x.data[i];

    sys->t += dt;

    free(k1); free(x_save);
    return true;
}

void etc_system_step_dynamic_trigger(ETCSystem* sys, double beta) {
    if (!sys) return;
    double x_norm_sq = etc_vector_dot(&sys->x, &sys->x);
    double e_norm_sq = etc_vector_dot(&sys->e, &sys->e);
    sys->eta_dot = -beta * sys->eta + sys->sigma * x_norm_sq - e_norm_sq;
    sys->eta += sys->dt * sys->eta_dot;
    if (sys->eta < 0.0) sys->eta = 0.0; /* Invariant: η(t) ≥ 0 */
}

void etc_system_simulate(ETCSystem* sys, double t_max, double dt) {
    if (!sys || t_max <= sys->t || dt <= 0.0) return;
    sys->dt = dt;
    sys->t_max = t_max;

    /* Record initial event (k=0) */
    etc_system_trigger_event(sys);

    while (sys->t < t_max) {
        if (!etc_system_step_rk4(sys)) break;

        /* Update dynamic trigger if applicable */
        if (sys->trigger_type == ETC_TRIGGER_DYNAMIC)
            etc_system_step_dynamic_trigger(sys, 1.0);

        /* Check trigger */
        if (etc_system_should_trigger(sys))
            etc_system_trigger_event(sys);

        /* Zeno check */
        if (sys->event_count > 10000) break;
    }
}

void etc_system_print(const ETCSystem* sys) {
    if (!sys) { printf("ETCSystem: NULL\n"); return; }
    printf("=== Event-Triggered Control System ===\n");
    printf("States: %d, Inputs: %d\n", sys->n_states, sys->n_inputs);
    printf("Trigger: type=%d, sigma=%.6f, epsilon=%.6f\n",
           sys->trigger_type, sys->sigma, sys->epsilon);
    printf("Time: t=%.4f, last_event=%.4f, dt=%.4f\n",
           sys->t, sys->t_last_event, sys->dt);
    printf("Events: %d, min IET=%.6f\n", sys->event_count, sys->min_iet_observed);
    printf("State: [");
    for (int i = 0; i < sys->n_states && i < 6; i++)
        printf("%.4f%s", sys->x.data[i], i < sys->n_states-1 ? ", " : "");
    if (sys->n_states > 6) printf(" ...");
    printf("]\n");
    printf("Error norm: %.6f, State norm: %.6f\n",
           etc_system_error_norm(sys), etc_system_state_norm(sys));
    printf("Regime: %d, Zeno bound: %.6f\n", sys->regime, sys->zeno_bound);
    printf("History: %d events, min_iet=%.6f, max_iet=%.6f, avg_iet=%.6f\n",
           sys->history.n_events, sys->history.min_iet,
           sys->history.max_iet, sys->history.avg_iet);
}

/* ============================================================================
 * Event History
 * ============================================================================ */

void etc_history_init(ETCEventHistory* hist, int capacity) {
    if (!hist) return;
    hist->capacity = capacity > 0 ? capacity : ETC_HISTORY_DEFAULT_CAP;
    hist->events = (ETCEvent*)calloc((size_t)hist->capacity, sizeof(ETCEvent));
    hist->n_events = 0;
    hist->min_iet = 1e100;
    hist->max_iet = 0.0;
    hist->avg_iet = 0.0;
    hist->zeno_detected = false;
}

void etc_history_append(ETCEventHistory* hist, const ETCEvent* event) {
    if (!hist || !event) return;
    /* Expand if needed */
    if (hist->n_events >= hist->capacity) {
        int new_cap = hist->capacity * 2;
        ETCEvent* new_ev = (ETCEvent*)realloc(hist->events,
                                              (size_t)new_cap * sizeof(ETCEvent));
        if (!new_ev) return;
        hist->events = new_ev;
        hist->capacity = new_cap;
    }
    /* Deep copy the event (vectors are owned by history after this) */
    ETCEvent* dst = &hist->events[hist->n_events];
    dst->time = event->time;
    dst->event_index = event->event_index;
    dst->trigger_value = event->trigger_value;
    dst->inter_event_time = event->inter_event_time;

    /* Deep copy vectors */
    int dim_state = event->state.dim;
    int dim_ctrl = event->control.dim;
    int dim_err = event->error_at_event.dim;
    dst->state = etc_vector_create(dim_state);
    dst->control = etc_vector_create(dim_ctrl);
    dst->error_at_event = etc_vector_create(dim_err);
    if (dst->state.data && event->state.data)
        for (int i = 0; i < dim_state; i++) dst->state.data[i] = event->state.data[i];
    if (dst->control.data && event->control.data)
        for (int i = 0; i < dim_ctrl; i++) dst->control.data[i] = event->control.data[i];
    if (dst->error_at_event.data && event->error_at_event.data)
        for (int i = 0; i < dim_err; i++) dst->error_at_event.data[i] = event->error_at_event.data[i];

    hist->n_events++;
}

void etc_history_compute_stats(ETCEventHistory* hist) {
    if (!hist || hist->n_events < 2) {
        if (hist) { hist->min_iet = 0.0; hist->max_iet = 0.0; hist->avg_iet = 0.0; }
        return;
    }
    double sum = 0.0;
    hist->min_iet = 1e100;
    hist->max_iet = 0.0;
    for (int i = 1; i < hist->n_events; i++) {
        double iet = hist->events[i].inter_event_time;
        if (iet < hist->min_iet) hist->min_iet = iet;
        if (iet > hist->max_iet) hist->max_iet = iet;
        sum += iet;
    }
    hist->avg_iet = sum / (double)(hist->n_events - 1);
}

void etc_history_free(ETCEventHistory* hist) {
    if (!hist) return;
    for (int i = 0; i < hist->n_events; i++) {
        etc_vector_free(&hist->events[i].state);
        etc_vector_free(&hist->events[i].control);
        etc_vector_free(&hist->events[i].error_at_event);
    }
    free(hist->events);
    hist->events = NULL;
    hist->n_events = 0;
    hist->capacity = 0;
}

bool etc_history_check_zeno(const ETCEventHistory* hist, double threshold) {
    if (!hist || hist->n_events < 10) return false;
    /* Check if event frequency exceeds threshold */
    if (hist->n_events >= 2) {
        double duration = hist->events[hist->n_events - 1].time
                        - hist->events[0].time;
        if (duration > 0.0) {
            double freq = (double)hist->n_events / duration;
            if (freq > threshold) return true;
        }
    }
    /* Check if minimum IET is below machine precision */
    if (hist->min_iet < 1e-12) return true;
    return false;
}
