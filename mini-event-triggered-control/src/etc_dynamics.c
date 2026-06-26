#include "etc_dynamics.h"
#include "etc_trigger.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Dynamics Kernel
 * ============================================================================ */

void etc_dynamics_closed_loop(const ETCSystem* sys, double* dx) {
    if (!sys || !dx) return;
    int n = sys->n_states;
    /* ẋ = Acl x − BK e */
    /* First compute Acl * x */
    ETCVector Aclx = etc_vector_create(n);
    etc_matrix_vec_mul(&sys->Acl, &sys->x, &Aclx);

    /* Compute BK * e */
    ETCMatrix BK = etc_matrix_create(n, n);
    etc_matrix_mul(&sys->B, &sys->K, &BK);
    ETCVector BKe = etc_vector_create(n);
    etc_matrix_vec_mul(&BK, &sys->e, &BKe);
    etc_matrix_free(&BK);

    for (int i = 0; i < n; i++)
        dx[i] = Aclx.data[i] - BKe.data[i];

    etc_vector_free(&Aclx);
    etc_vector_free(&BKe);
}

void etc_dynamics_open_loop(const ETCSystem* sys, double* dx) {
    if (!sys || !dx) return;
    etc_system_compute_derivative(sys, dx);
}

/* ============================================================================
 * Matrix Exponential (Scaling-and-Squaring + Taylor)
 *
 * e^{At} ≈ (e^{At / 2^s})^{2^s}
 * where e^{At / 2^s} is computed via truncated Taylor series.
 * ============================================================================ */

void etc_matrix_exponential(const ETCMatrix* A, double t, int n,
                             ETCMatrix* expAt) {
    if (!A || !expAt || n <= 0 || !A->data || !expAt->data) return;

    int n2 = n * n;

    /* Scaling: choose s such that ||A t / 2^s|| < 1 */
    double norm_A = etc_matrix_norm_2(A);
    double scaled_norm = norm_A * fabs(t);
    int s = 0;
    while (scaled_norm > 1.0 && s < 20) {
        scaled_norm /= 2.0;
        s++;
    }

    /* Scale A: A_scaled = A * t / 2^s */
    double* A_scaled = (double*)malloc((size_t)n2 * sizeof(double));
    double scale = t / (double)(1 << s);
    for (int i = 0; i < n2; i++) A_scaled[i] = A->data[i] * scale;

    /* Taylor series: e^{A_scaled} ≈ I + A_scaled + A_scaled^2/2! + ... */
    double* Term = (double*)malloc((size_t)n2 * sizeof(double));
    double* Sum = (double*)malloc((size_t)n2 * sizeof(double));
    double* Temp = (double*)malloc((size_t)n2 * sizeof(double));
    if (!A_scaled || !Term || !Sum || !Temp) {
        free(A_scaled); free(Term); free(Sum); free(Temp);
        return;
    }

    /* Initialize: Sum = I, Term = I */
    for (int i = 0; i < n2; i++) {
        Sum[i] = (i % (n + 1) == 0) ? 1.0 : 0.0;
        Term[i] = Sum[i];
    }

    int degree = 15;
    for (int k = 1; k <= degree; k++) {
        /* Term = Term * A_scaled / k */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double sval = 0.0;
                for (int p = 0; p < n; p++)
                    sval += Term[i * n + p] * A_scaled[p * n + j];
                Temp[i * n + j] = sval / (double)k;
            }
        }
        /* Add to Sum and update Term */
        for (int i = 0; i < n2; i++) {
            Term[i] = Temp[i];
            Sum[i] += Term[i];
        }
    }

    /* Squaring: Sum = (Sum)^{2^s} */
    for (int step = 0; step < s; step++) {
        /* Temp = Sum */
        for (int i = 0; i < n2; i++) Temp[i] = Sum[i];
        /* Sum = Temp * Temp */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double sval = 0.0;
                for (int p = 0; p < n; p++)
                    sval += Temp[i * n + p] * Temp[p * n + j];
                Sum[i * n + j] = sval;
            }
        }
    }

    /* Copy result */
    for (int i = 0; i < n2; i++) expAt->data[i] = Sum[i];

    free(A_scaled); free(Term); free(Sum); free(Temp);
}

/* ============================================================================
 * Flow Map: Φ_t(x₀, u) = e^{At} x₀ + ∫₀ᵗ e^{A(t−s)} B ds · u
 *
 * For constant u (zero-order hold), the integral term simplifies:
 *   ∫₀ᵗ e^{A(t−s)} B ds = A^{-1} (e^{At} − I) B   (if A invertible)
 *
 * For general A (possibly singular), we compute the augmented system:
 *   d/dt [x; 1] = [A  Bu; 0  0] [x; 1]
 *   with initial condition [x₀; 1].
 * ============================================================================ */

void etc_flow_map(const ETCMatrix* A, const ETCMatrix* B,
                   const ETCVector* x0, const ETCVector* u_const,
                   double t, int n, int m, ETCVector* x_t) {
    if (!A || !B || !x0 || !u_const || !x_t) return;
    if (n <= 0 || m <= 0) return;

    /* Compute Bu */
    ETCVector Bu = etc_vector_create(n);
    etc_matrix_vec_mul(B, u_const, &Bu);

    /* For small n: use augmented system approach.
     * Augmented matrix M = [A  Bu; 0  0] of size (n+1)×(n+1)
     * e^{M t} * [x0; 1] gives [x(t); 1] */
    int n_aug = n + 1;
    ETCMatrix M = etc_matrix_create(n_aug, n_aug);

    /* M = [A  Bu] */
    /*     [0   0] */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            M.data[i * n_aug + j] = A->data[i * n + j];
        M.data[i * n_aug + n] = Bu.data[i]; /* Bu in last column */
    }
    /* Last row is already zero (calloc) */

    ETCMatrix expM = etc_matrix_create(n_aug, n_aug);
    etc_matrix_exponential(&M, t, n_aug, &expM);

    /* x(t) = expM[0:n, 0:n] * x0 + expM[0:n, n] * 1 */
    for (int i = 0; i < n; i++) {
        x_t->data[i] = 0.0;
        for (int j = 0; j < n; j++)
            x_t->data[i] += expM.data[i * n_aug + j] * x0->data[j];
        x_t->data[i] += expM.data[i * n_aug + n]; /* times 1 */
    }

    etc_vector_free(&Bu);
    etc_matrix_free(&M);
    etc_matrix_free(&expM);
}

/* ============================================================================
 * Simulation Engine
 * ============================================================================ */

void etc_sim_config_init(ETCSimConfig* cfg) {
    if (!cfg) return;
    cfg->t_start = 0.0;
    cfg->t_end = 10.0;
    cfg->dt = 0.001;
    cfg->dt_max = 0.01;
    cfg->dt_min = 1e-6;
    cfg->max_events = 100000;
    cfg->max_steps = 10000000;
    cfg->detect_zeno = true;
    cfg->zeno_threshold = 1e6;
    cfg->event_tolerance = 1e-8;
}

void etc_simulate_full(const ETCSystem* sys_const, const ETCSimConfig* cfg,
                        ETCSimResult* result) {
    if (!sys_const || !cfg || !result) return;

    /* Work on a mutable copy */
    ETCSystem* sys = (ETCSystem*)malloc(sizeof(ETCSystem));
    if (!sys) return;
    memcpy(sys, sys_const, sizeof(ETCSystem));

    memset(result, 0, sizeof(ETCSimResult));

    double t = cfg->t_start;
    int step_count = 0;
    int event_count = 0;
    double* iet_buffer = NULL;
    int iet_cap = 0;
    int iet_count = 0;
    double sum_error = 0.0;
    double sum_state = 0.0;
    int sample_count = 0;
    double settling_time = cfg->t_end;
    double x0_norm = etc_vector_norm(&sys->x);
    bool settled = false;

    /* Record initial event */
    etc_system_trigger_event(sys);
    event_count = 1;

    while (t < cfg->t_end && step_count < cfg->max_steps
           && event_count < cfg->max_events) {
        /* Integrate one step */
        bool event_fired = false;
        if (!etc_step_integrate(sys, cfg->dt, &event_fired)) break;
        t = sys->t;
        step_count++;

        if (event_fired) {
            event_count++;
            /* Record IET */
            double iet = sys->history.events[sys->history.n_events - 1].inter_event_time;
            if (iet_count >= iet_cap) {
                int new_cap = iet_cap == 0 ? 256 : iet_cap * 2;
                iet_buffer = (double*)realloc(iet_buffer,
                                              (size_t)new_cap * sizeof(double));
                iet_cap = new_cap;
            }
            if (iet_buffer) iet_buffer[iet_count++] = iet;
        }

        /* Accumulate statistics */
        sum_error += etc_system_error_norm(sys);
        sum_state += etc_system_state_norm(sys);
        sample_count++;

        /* Settling time check */
        if (!settled) {
            double xn = etc_vector_norm(&sys->x);
            if (xn < 0.02 * x0_norm + 1e-6) {
                settling_time = t;
                settled = true;
            }
        }
    }

    /* Fill results */
    result->t_final = t;
    result->n_steps = step_count;
    result->n_events = event_count;
    result->settling_time = settling_time;
    result->converged = settled;
    result->total_control_updates = (double)event_count;

    if (iet_count > 0 && iet_buffer) {
        double min_iet = iet_buffer[0], max_iet = iet_buffer[0], sum_iet = 0.0;
        for (int i = 0; i < iet_count; i++) {
            double v = iet_buffer[i];
            if (v < min_iet) min_iet = v;
            if (v > max_iet) max_iet = v;
            sum_iet += v;
        }
        result->min_iet = min_iet;
        result->max_iet = max_iet;
        result->avg_iet = sum_iet / (double)iet_count;

        /* Standard deviation */
        double mean = result->avg_iet;
        double var = 0.0;
        for (int i = 0; i < iet_count; i++) {
            double d = iet_buffer[i] - mean;
            var += d * d;
        }
        result->std_iet = sqrt(var / (double)iet_count);
    }

    if (sample_count > 0) {
        result->avg_error_norm = sum_error / (double)sample_count;
        result->avg_state_norm = sum_state / (double)sample_count;
    }

    result->zeno_free = !etc_history_check_zeno(&sys->history, cfg->zeno_threshold);

    free(iet_buffer);
    /* Don't free sys's internal data — it's shared with sys_const */
    free(sys);
}

bool etc_step_integrate(ETCSystem* sys, double dt, bool* event_fired) {
    if (!sys) return false;
    sys->dt = dt;

    /* Evaluate trigger before step */
    double gamma_before = etc_system_eval_trigger(sys);
    bool was_near = (gamma_before > -1e-6);

    /* Integrate */
    if (!etc_system_step_rk4(sys)) return false;

    /* Evaluate trigger after step */
    double gamma_after = etc_system_eval_trigger(sys);
    bool is_over = (gamma_after >= 0.0);

    if (event_fired) *event_fired = false;

    if (is_over) {
        /* Event triggered: reset sampled state */
        etc_system_trigger_event(sys);
        if (event_fired) *event_fired = true;
    } else if (was_near && gamma_after < gamma_before) {
        /* Approaching trigger from above: could interpolate for better accuracy */
        /* For now, standard step is sufficient */
    }

    return true;
}

bool etc_detect_event(ETCSystem* sys, double dt, double* event_time) {
    if (!sys || !event_time) return false;

    double t0 = sys->t;
    double gamma0 = etc_system_eval_trigger(sys);

    /* Save state */
    int n = sys->n_states;
    double* x_save = (double*)malloc((size_t)n * sizeof(double));
    double eta_save = sys->eta;
    for (int i = 0; i < n; i++) x_save[i] = sys->x.data[i];

    /* Take a trial step */
    sys->dt = dt;
    if (!etc_system_step_rk4(sys)) {
        free(x_save);
        return false;
    }

    double gamma1 = etc_system_eval_trigger(sys);
    bool crossing = (gamma0 < 0.0 && gamma1 >= 0.0);

    /* Restore state */
    for (int i = 0; i < n; i++) sys->x.data[i] = x_save[i];
    sys->eta = eta_save;
    sys->t = t0;
    /* Restore error */
    for (int i = 0; i < n; i++)
        sys->e.data[i] = sys->x_hat.data[i] - sys->x.data[i];

    if (crossing) {
        /* Linear interpolation for crossing time */
        double alpha = -gamma0 / (gamma1 - gamma0 + 1e-15);
        *event_time = t0 + alpha * dt;
        free(x_save);
        return true;
    }

    free(x_save);
    return false;
}

/* ============================================================================
 * Inter-Event Time Analysis
 * ============================================================================ */

double etc_inter_event_time_bound(const ETCSystem* sys) {
    if (!sys) return 0.0;
    int n = sys->n_states;
    ETCMatrix BK = etc_matrix_create(n, n);
    etc_matrix_mul(&sys->B, &sys->K, &BK);
    double bound = etc_compute_iet_lower_bound(&sys->Acl, &BK, sys->sigma);
    etc_matrix_free(&BK);
    return bound;
}

double etc_iet_estimate(const ETCSystem* sys, const ETCVector* xk) {
    if (!sys || !xk) return 0.0;
    int n = sys->n_states;

    /* The error dynamics between events (ė = −Acl x − BK e) indicate
     * that the time to violation depends on |Acl x_k| and the threshold.
     *
     * A simple estimate: τ ≈ σ|x_k| / |Acl x_k|
     *   (time for error to grow to σ|x| under worst-case dynamics)
     */
    ETCVector Acl_xk = etc_vector_create(n);
    etc_matrix_vec_mul(&sys->Acl, xk, &Acl_xk);
    double norm_Aclxk = etc_vector_norm(&Acl_xk);
    double norm_xk = etc_vector_norm(xk);
    etc_vector_free(&Acl_xk);

    if (norm_Aclxk < 1e-12) return sys->t_max - sys->t; /* No more events needed */

    double tau = sys->sigma * norm_xk / norm_Aclxk;
    if (tau < 1e-8) tau = 1e-8;
    return tau;
}

double etc_average_iet(const ETCEventHistory* hist) {
    if (!hist || hist->n_events < 2) return 0.0;
    double sum = 0.0;
    for (int i = 1; i < hist->n_events; i++)
        sum += hist->events[i].inter_event_time;
    return sum / (double)(hist->n_events - 1);
}

double etc_communication_ratio(int n_events, double t_total, double dt) {
    if (t_total <= 0.0 || dt <= 0.0) return 1.0;
    double periodic_samples = t_total / dt;
    if (periodic_samples < 1.0) return 1.0;
    return (double)n_events / periodic_samples;
}

void etc_sim_result_print(const ETCSimResult* result) {
    if (!result) return;
    printf("=== ETC Simulation Results ===\n");
    printf("Final time:    %.4f s\n", result->t_final);
    printf("Steps:         %d\n", result->n_steps);
    printf("Events:        %d\n", result->n_events);
    printf("Min IET:       %.6f s\n", result->min_iet);
    printf("Max IET:       %.6f s\n", result->max_iet);
    printf("Avg IET:       %.6f s\n", result->avg_iet);
    printf("Std IET:       %.6f s\n", result->std_iet);
    printf("Avg |e|:       %.6f\n", result->avg_error_norm);
    printf("Avg |x|:       %.6f\n", result->avg_state_norm);
    printf("Settling time: %.4f s\n", result->settling_time);
    printf("Converged:     %s\n", result->converged ? "YES" : "NO");
    printf("Zeno-free:     %s\n", result->zeno_free ? "YES" : "NO");
    printf("Control updates: %.0f\n", result->total_control_updates);
}

/* ============================================================================
 * Adaptive Step-Size Control
 * ============================================================================ */

bool etc_adaptive_step(const ETCSystem* sys, double dt_current,
                        double* dt_next) {
    if (!sys || !dt_next) return false;

    double gamma = etc_system_eval_trigger(sys);
    double margin = -gamma; /* Positive = safe, negative = already triggered */

    /* PID-like step size adjustment based on trigger margin */
    double safety_factor = 0.9;
    double dt_min = 1e-6;
    double dt_max = 0.1;

    if (margin > 1.0) {
        /* Far from trigger: increase step */
        *dt_next = dt_current * 1.5;
    } else if (margin > 0.1) {
        /* Moderate margin: slight increase */
        *dt_next = dt_current * 1.1;
    } else if (margin > 0.01) {
        /* Near trigger: maintain step */
        *dt_next = dt_current;
    } else if (margin > 0.0) {
        /* Very near trigger: reduce step */
        *dt_next = dt_current * 0.5;
    } else {
        /* Already triggered or past: minimum step */
        *dt_next = dt_min;
    }

    if (*dt_next > dt_max) *dt_next = dt_max;
    if (*dt_next < dt_min) *dt_next = dt_min;
    *dt_next *= safety_factor;

    return true;
}
