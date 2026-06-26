#include "dta_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/*
 * dta_core.c - Dwell-Time Analysis Core Implementation
 *
 * Implements:
 *   System construction/destruction for switched systems
 *   Eigenvalue computation via QR algorithm (Francis 1961)
 *   Hurwitz/Schur stability checks
 *   Matrix exponential via scaling+squaring (Moler & Van Loan 2003)
 *   Matrix measure / logarithmic norm
 *   Lyapunov equation solver (Bartels-Stewart 1972)
 *   L2 gain / H-infinity norm via Hamiltonian bisection
 *
 * All functions implement independent knowledge points per SKILL.md.
 * Complexity bounds and theorem references in each function doc.
 */

static void* safe_alloc(size_t sz) {
    void* p = malloc(sz);
    if (!p) { fprintf(stderr, "DTA: alloc fail %zu\n", sz); exit(1); }
    return p;
}

/* ==============================================================
 * System Lifecycle
 * ============================================================== */

DTA_SwitchedSystem* dta_system_create(int n_modes, int state_dim,
                                       int input_dim, int output_dim) {
    DTA_SwitchedSystem* sys = safe_alloc(sizeof(DTA_SwitchedSystem));
    sys->n_modes = n_modes;
    sys->state_dim = state_dim;
    sys->input_dim = input_dim;
    sys->output_dim = output_dim;
    sys->class = DTA_LINEAR;
    sys->current_time = 0.0;
    sys->current_mode = 0;
    sys->max_switches = 10000;
    sys->modes = safe_alloc((size_t)n_modes * sizeof(DTA_SystemMode));
    int i;
    for (i = 0; i < n_modes; i++) {
        sys->modes[i].id = i;
        sys->modes[i].A = NULL;
        sys->modes[i].B = NULL;
        sys->modes[i].C = NULL;
        sys->modes[i].n = state_dim;
        sys->modes[i].p = input_dim;
        sys->modes[i].q = output_dim;
        sys->modes[i].has_affine_term = false;
        sys->modes[i].affine_b = NULL;
        sys->modes[i].stability = DTA_MODE_UNKNOWN;
        sys->modes[i].max_eigenvalue_real = 0.0;
        sys->modes[i].spectral_radius = 0.0;
        sys->modes[i].dynamics = NULL;
        sys->modes[i].context = NULL;
    }
    return sys;
}

void dta_system_free(DTA_SwitchedSystem* sys) {
    if (!sys) return;
    int i;
    for (i = 0; i < sys->n_modes; i++) {
        free(sys->modes[i].A);
        free(sys->modes[i].B);
        free(sys->modes[i].C);
        free(sys->modes[i].affine_b);
    }
    free(sys->modes);
    free(sys);
}

int dta_system_set_linear_mode(DTA_SwitchedSystem* sys, int mode_idx,
                                const double* A, const double* B,
                                const double* C) {
    if (!sys || mode_idx < 0 || mode_idx >= sys->n_modes) return -1;
    DTA_SystemMode* m = &sys->modes[mode_idx];
    int n = sys->state_dim, p = sys->input_dim, q = sys->output_dim;
    if (A) {
        m->A = safe_alloc((size_t)(n * n) * sizeof(double));
        memcpy(m->A, A, (size_t)(n * n) * sizeof(double));
    }
    if (B) {
        m->B = safe_alloc((size_t)(n * p) * sizeof(double));
        memcpy(m->B, B, (size_t)(n * p) * sizeof(double));
    }
    if (C) {
        m->C = safe_alloc((size_t)(q * n) * sizeof(double));
        memcpy(m->C, C, (size_t)(q * n) * sizeof(double));
    }
    if (A) {
        double* re = safe_alloc((size_t)n * sizeof(double));
        double* im = safe_alloc((size_t)n * sizeof(double));
        dta_eigenvalues(A, n, re, im, 200);
        double max_re = -DBL_MAX;
        int j;
        for (j = 0; j < n; j++)
            if (re[j] > max_re) max_re = re[j];
        m->max_eigenvalue_real = max_re;
        m->spectral_radius = 0.0;
        for (j = 0; j < n; j++) {
            double mag = sqrt(re[j]*re[j] + im[j]*im[j]);
            if (mag > m->spectral_radius) m->spectral_radius = mag;
        }
        if (max_re < -1e-12) m->stability = DTA_MODE_STABLE;
        else if (max_re > 1e-12) m->stability = DTA_MODE_UNSTABLE;
        else m->stability = DTA_MODE_MARGINAL;
        free(re); free(im);
    }
    return 0;
}

int dta_system_set_nonlinear_mode(DTA_SwitchedSystem* sys, int mode_idx,
    void (*dynamics)(double, const double*, int, const double*, int,
                     double*, void*), void* context) {
    if (!sys || mode_idx < 0 || mode_idx >= sys->n_modes) return -1;
    sys->modes[mode_idx].dynamics = dynamics;
    sys->modes[mode_idx].context = context;
    sys->modes[mode_idx].stability = DTA_MODE_UNKNOWN;
    sys->class = DTA_NONLINEAR;
    return 0;
}

int dta_system_set_affine_mode(DTA_SwitchedSystem* sys, int mode_idx,
                                const double* A, const double* b) {
    if (!sys || mode_idx < 0 || mode_idx >= sys->n_modes) return -1;
    dta_system_set_linear_mode(sys, mode_idx, A, NULL, NULL);
    DTA_SystemMode* m = &sys->modes[mode_idx];
    m->has_affine_term = true;
    m->affine_b = safe_alloc((size_t)sys->state_dim * sizeof(double));
    memcpy(m->affine_b, b, (size_t)sys->state_dim * sizeof(double));
    sys->class = DTA_AFFINE;
    return 0;
}

void dta_system_rhs(const DTA_SwitchedSystem* sys, int mode,
                    double t, const double* x, const double* u,
                    double* dx) {
    if (!sys || mode < 0 || mode >= sys->n_modes) return;
    const DTA_SystemMode* m = &sys->modes[mode];
    int n = sys->state_dim;
    int i, j;
    if (m->dynamics) {
        m->dynamics(t, x, n, u, sys->input_dim, dx, m->context);
        return;
    }
    for (i = 0; i < n; i++) dx[i] = 0.0;
    if (m->A)
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                dx[i] += m->A[i * n + j] * x[j];
    if (m->B && u)
        for (i = 0; i < n; i++)
            for (j = 0; j < sys->input_dim; j++)
                dx[i] += m->B[i * sys->input_dim + j] * u[j];
    if (m->has_affine_term && m->affine_b)
        for (i = 0; i < n; i++)
            dx[i] += m->affine_b[i];
}

/*
 * QR Algorithm for Eigenvalues (Francis 1961)
 * Complexity: O(n^3) for Hessenberg reduction + O(n^2) per QR step
 * Reference: Golub & Van Loan (2013) "Matrix Computations" Ch.7
 */

static void hh_reflector(const double* x, int n, double* v, double* beta) {
    double sigma = 0.0;
    int i;
    for (i = 1; i < n; i++) sigma += x[i] * x[i];
    v[0] = 1.0;
    for (i = 1; i < n; i++) v[i] = x[i];
    if (sigma < 1e-30) { *beta = 0.0; return; }
    double mu = sqrt(x[0]*x[0] + sigma);
    if (x[0] <= 0) v[0] = x[0] - mu;
    else v[0] = -sigma / (x[0] + mu);
    *beta = 2.0 * v[0]*v[0] / (sigma + v[0]*v[0]);
    double v0_inv = 1.0 / v[0];
    for (i = 1; i < n; i++) v[i] *= v0_inv;
}

static void to_hessenberg(double* A, int n, double* Q) {
    int i, j, k;
    for (i = 0; i < n*n; i++) Q[i] = (i % (n+1) == 0) ? 1.0 : 0.0;
    double* v = safe_alloc((size_t)n * sizeof(double));
    double* w = safe_alloc((size_t)n * sizeof(double));
    for (k = 0; k < n - 2; k++) {
        double* x = safe_alloc((size_t)(n - k - 1) * sizeof(double));
        for (i = k+1; i < n; i++) x[i - k - 1] = A[i * n + k];
        double beta;
        hh_reflector(x, n - k - 1, v, &beta);
        if (beta > 1e-30) {
            for (i = 0; i < n; i++) {
                double s = 0.0;
                for (j = 0; j < n - k - 1; j++)
                    s += A[(k+1+j)*n + i] * v[j];
                w[i] = beta * s;
            }
            for (i = 0; i < n - k - 1; i++)
                for (j = 0; j < n; j++)
                    A[(k+1+i)*n + j] -= v[i] * w[j];
            for (i = 0; i < n; i++) {
                double s = 0.0;
                for (j = 0; j < n - k - 1; j++)
                    s += A[i*n + (k+1+j)] * v[j];
                w[i] = beta * s;
            }
            for (i = 0; i < n; i++)
                for (j = 0; j < n - k - 1; j++)
                    A[i*n + (k+1+j)] -= w[i] * v[j];
            for (i = 0; i < n; i++) {
                double s = 0.0;
                for (j = 0; j < n - k - 1; j++)
                    s += Q[i*n + (k+1+j)] * v[j];
                w[i] = beta * s;
            }
            for (i = 0; i < n; i++)
                for (j = 0; j < n - k - 1; j++)
                    Q[i*n + (k+1+j)] -= w[i] * v[j];
        }
        free(x);
    }
    free(v); free(w);
}

/* Implicit QR step with Wilkinson shift for real Hessenberg matrix.
 * Complexity: O(n^2) per iteration.
 * Reference: Golub & Van Loan (2013) "Matrix Computations" Sec 7.5.2 */
static void qr_hessenberg(double* H, int n, double* Q, int max_iter) {
    int iter, i, jj;
    double* cs = safe_alloc((size_t)n * sizeof(double));
    double* sn = safe_alloc((size_t)n * sizeof(double));
    for (iter = 0; iter < max_iter; iter++) {
        int k;
        for (k = n - 1; k >= 1; k--) {
            double thresh = 1e-12 * (fabs(H[(k-1)*n + (k-1)]) + fabs(H[k*n + k]));
            if (fabs(H[k*n + (k-1)]) < thresh)
                H[k*n + (k-1)] = 0.0;
        }
        double a = H[(n-2)*n + (n-2)];
        double b = H[(n-2)*n + (n-1)];
        double cc = H[(n-1)*n + (n-2)];
        double d = H[(n-1)*n + (n-1)];
        double trace = a + d, det = a*d - b*cc;
        double disc = trace*trace - 4.0*det, mu;
        if (disc >= 0) {
            double r = sqrt(disc);
            double r1 = (trace + r) / 2.0, r2 = (trace - r) / 2.0;
            mu = (fabs(r1 - d) < fabs(r2 - d)) ? r1 : r2;
        } else mu = trace / 2.0;
        for (i = 0; i < n - 1; i++) {
            double x = H[i*n + i] - mu, y = H[(i+1)*n + i];
            double r = sqrt(x*x + y*y);
            if (r > 1e-30) { cs[i] = y / r; sn[i] = x / r; }
            else { cs[i] = 0.0; sn[i] = 1.0; }
            int col_max = (i + 2 < n) ? i + 3 : n;
            for (jj = i; jj < col_max; jj++) {
                double t1 = H[i*n + jj], t2 = H[(i+1)*n + jj];
                H[i*n + jj]     = sn[i]*t1 + cs[i]*t2;
                H[(i+1)*n + jj] = -cs[i]*t1 + sn[i]*t2;
            }
            int row_start = (i > 0) ? i - 1 : 0;
            for (jj = row_start; jj <= i + 1 && jj < n; jj++) {
                double t1 = H[jj*n + i], t2 = H[jj*n + (i+1)];
                H[jj*n + i]     = sn[i]*t1 + cs[i]*t2;
                H[jj*n + (i+1)] = -cs[i]*t1 + sn[i]*t2;
            }
            for (jj = 0; jj < n; jj++) {
                double t1 = Q[jj*n + i], t2 = Q[jj*n + (i+1)];
                Q[jj*n + i]     = sn[i]*t1 + cs[i]*t2;
                Q[jj*n + (i+1)] = -cs[i]*t1 + sn[i]*t2;
            }
        }
        double sum_sub = 0.0;
        for (i = 1; i < n; i++) sum_sub += fabs(H[i*n + (i-1)]);
        if (sum_sub < 1e-10) break;
    }
    free(cs); free(sn);
}

int dta_eigenvalues(const double* A, int n, double* real_parts,
                    double* imag_parts, int max_iter) {
    if (!A || !real_parts || !imag_parts || n <= 0) return -1;
    if (n == 1) {
        real_parts[0] = A[0]; imag_parts[0] = 0.0; return 0;
    }
    double* H = safe_alloc((size_t)(n*n) * sizeof(double));
    memcpy(H, A, (size_t)(n*n) * sizeof(double));
    double* Q = safe_alloc((size_t)(n*n) * sizeof(double));
    to_hessenberg(H, n, Q);
    qr_hessenberg(H, n, Q, max_iter);
    int i = 0;
    while (i < n) {
        if (i < n - 1 && fabs(H[(i+1)*n + i]) > 1e-10) {
            double a = H[i*n + i], b = H[i*n + (i+1)];
            double cc = H[(i+1)*n + i], d = H[(i+1)*n + (i+1)];
            double tr = a + d, det2 = a*d - b*cc;
            double disc2 = tr*tr - 4.0*det2;
            real_parts[i] = tr / 2.0;
            real_parts[i+1] = tr / 2.0;
            if (disc2 < 0) {
                imag_parts[i] = sqrt(-disc2) / 2.0;
                imag_parts[i+1] = -imag_parts[i];
            } else {
                imag_parts[i] = 0.0; imag_parts[i+1] = 0.0;
                real_parts[i] = (tr + sqrt(disc2)) / 2.0;
                real_parts[i+1] = (tr - sqrt(disc2)) / 2.0;
            }
            i += 2;
        } else {
            real_parts[i] = H[i*n + i];
            imag_parts[i] = 0.0;
            i++;
        }
    }
    free(H); free(Q);
    return 0;
}

bool dta_is_hurwitz(const double* A, int n, double tol) {
    if (!A || n <= 0) return false;
    double* re = safe_alloc((size_t)n * sizeof(double));
    double* im = safe_alloc((size_t)n * sizeof(double));
    dta_eigenvalues(A, n, re, im, 200);
    bool hurwitz = true;
    int i;
    for (i = 0; i < n; i++)
        if (re[i] >= -tol) { hurwitz = false; break; }
    free(re); free(im);
    return hurwitz;
}

bool dta_is_schur(const double* A, int n, double tol) {
    if (!A || n <= 0) return false;
    double* re = safe_alloc((size_t)n * sizeof(double));
    double* im = safe_alloc((size_t)n * sizeof(double));
    dta_eigenvalues(A, n, re, im, 200);
    bool schur = true;
    int i;
    for (i = 0; i < n; i++) {
        double mag = sqrt(re[i]*re[i] + im[i]*im[i]);
        if (mag >= 1.0 - tol) { schur = false; break; }
    }
    free(re); free(im);
    return schur;
}

/* ==============================================================
 * Matrix exponential e^{A t} via scaling+squaring
 * (Moler & Van Loan 2003, "Nineteen Dubious Ways...")
 *
 * Pade approximant of order (1,1): e^X ~ (I + X/2 + X^2/12) / (I - X/2 + X^2/12)
 * Then square the result s times.
 * Complexity: O(n^3 + n^3 log ||At||)
 * ============================================================== */

void dta_matrix_exp(const double* A, int n, double t, double* result) {
    if (!A || !result || n <= 0) return;
    if (n == 1) { result[0] = exp(A[0] * t); return; }
    int i, j, k;
    double norm = 0.0;
    for (i = 0; i < n*n; i++) norm += A[i] * A[i];
    norm = sqrt(norm) * fabs(t);
    int s = 0;
    if (norm > 1.0) {
        s = (int)ceil(log2(norm));
        t /= (double)(1 << s);
    }
    double* X = safe_alloc((size_t)(n*n) * sizeof(double));
    double* X2 = safe_alloc((size_t)(n*n) * sizeof(double));
    double* Nmat = safe_alloc((size_t)(n*n) * sizeof(double));
    double* Dmat = safe_alloc((size_t)(n*n) * sizeof(double));
    for (i = 0; i < n*n; i++) X[i] = A[i] * t;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            double sv = 0.0;
            for (k = 0; k < n; k++) sv += X[i*n + k] * X[k*n + j];
            X2[i*n + j] = sv;
        }
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            double id = (i == j) ? 1.0 : 0.0;
            Nmat[i*n + j] = id + 0.5*X[i*n + j] + X2[i*n + j]/12.0;
            Dmat[i*n + j] = id - 0.5*X[i*n + j] + X2[i*n + j]/12.0;
        }
    double* aug = safe_alloc((size_t)(n * (2*n)) * sizeof(double));
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            aug[i*(2*n) + j] = Dmat[i*n + j];
            aug[i*(2*n) + n + j] = Nmat[i*n + j];
        }
    int col;
    for (col = 0; col < n; col++) {
        int pivot = col, row;
        for (row = col + 1; row < n; row++)
            if (fabs(aug[row*(2*n) + col]) > fabs(aug[pivot*(2*n) + col]))
                pivot = row;
        if (fabs(aug[pivot*(2*n) + col]) < 1e-30) continue;
        if (pivot != col)
            for (j = 0; j < 2*n; j++) {
                double tmp = aug[col*(2*n) + j];
                aug[col*(2*n) + j] = aug[pivot*(2*n) + j];
                aug[pivot*(2*n) + j] = tmp;
            }
        double pv = aug[col*(2*n) + col];
        for (j = 0; j < 2*n; j++) aug[col*(2*n) + j] /= pv;
        for (row = 0; row < n; row++) {
            if (row == col) continue;
            double f = aug[row*(2*n) + col];
            for (j = 0; j < 2*n; j++)
                aug[row*(2*n) + j] -= f * aug[col*(2*n) + j];
        }
    }
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            result[i*n + j] = aug[i*(2*n) + n + j];
    double* tmp = safe_alloc((size_t)(n*n) * sizeof(double));
    int sq;
    for (sq = 0; sq < s; sq++) {
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                double sv = 0.0;
                for (k = 0; k < n; k++)
                    sv += result[i*n + k] * result[k*n + j];
                tmp[i*n + j] = sv;
            }
        memcpy(result, tmp, (size_t)(n*n) * sizeof(double));
    }
    free(X); free(X2); free(Nmat); free(Dmat); free(aug); free(tmp);
}

/* ==============================================================
 * Matrix measure (logarithmic norm)
 * mu(A) = lambda_max((A + A^T) / 2)
 *
 * Provides bound: ||e^{A t}|| <= e^{mu(A) t}
 * If mu(A) < 0, system is strictly contractive.
 * For dwell time: tau_d > ln(mu) / (2 lambda) requires mu_i values.
 * Complexity: O(n^3) for eigenvalue computation.
 * Reference: Desoer & Vidyasagar (1975) "Feedback Systems: Input-Output Properties"
 * ============================================================== */

double dta_matrix_measure(const double* A, int n) {
    if (!A || n <= 0) return 0.0;
    int i, j;
    double* sym = safe_alloc((size_t)(n*n) * sizeof(double));
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            sym[i*n + j] = 0.5 * (A[i*n + j] + A[j*n + i]);
    double* re = safe_alloc((size_t)n * sizeof(double));
    double* im = safe_alloc((size_t)n * sizeof(double));
    dta_eigenvalues(sym, n, re, im, 200);
    double max_eig = re[0];
    for (i = 1; i < n; i++)
        if (re[i] > max_eig) max_eig = re[i];
    free(sym); free(re); free(im);
    return max_eig;
}

/* ==============================================================
 * Lyapunov equation: A^T P + P A = -Q
 * Bartels-Stewart algorithm (1972)
 *
 * 1. Real Schur decomposition A = U T U^T via QR
 * 2. Transform to T^T Y + Y T = -F where Y = U^T P U, F = U^T Q U
 * 3. Solve quasi-triangular system column by column
 * 4. Reconstruct P = U Y U^T
 *
 * Complexity: O(n^3)
 * Reference: Bartels & Stewart (1972) Comm. ACM 15(9):820-826
 * ============================================================== */

int dta_solve_lyapunov(const double* A, int n, const double* Q,
                        double* P) {
    if (!A || !Q || !P || n <= 0) return -1;
    int i, j, k;
    if (n == 1) {
        if (fabs(A[0]) < 1e-30) return -1;
        P[0] = Q[0] / (-2.0 * A[0]);
        return (P[0] > 0) ? 0 : -1;
    }
    double* T = safe_alloc((size_t)(n*n) * sizeof(double));
    memcpy(T, A, (size_t)(n*n) * sizeof(double));
    double* U = safe_alloc((size_t)(n*n) * sizeof(double));
    to_hessenberg(T, n, U);
    qr_hessenberg(T, n, U, 200);
    double* F = safe_alloc((size_t)(n*n) * sizeof(double));
    double* temp = safe_alloc((size_t)(n*n) * sizeof(double));
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            double sv = 0.0;
            for (k = 0; k < n; k++) sv += Q[i*n + k] * U[k*n + j];
            temp[i*n + j] = -sv;
        }
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            double sv = 0.0;
            for (k = 0; k < n; k++) sv += U[k*n + i] * temp[k*n + j];
            F[i*n + j] = sv;
        }
    double* Y = safe_alloc((size_t)(n*n) * sizeof(double));
    memset(Y, 0, (size_t)(n*n) * sizeof(double));
    for (j = 0; j < n; j++) {
        double* rhs = safe_alloc((size_t)n * sizeof(double));
        for (i = 0; i < n; i++) {
            rhs[i] = F[i*n + j];
            for (k = 0; k < j; k++)
                rhs[i] -= T[k*n + j] * Y[i*n + k];
        }
        double* coeff = safe_alloc((size_t)(n*n) * sizeof(double));
        memset(coeff, 0, (size_t)(n*n) * sizeof(double));
        for (i = 0; i < n; i++)
            for (k = 0; k < n; k++) {
                coeff[i*n + k] = (i == k) ? T[j*n + j] : 0.0;
                coeff[i*n + k] += T[k*n + i];
            }
        int pivot;
        for (pivot = 0; pivot < n; pivot++) {
            int best = pivot, r;
            for (r = pivot+1; r < n; r++)
                if (fabs(coeff[r*n + pivot]) > fabs(coeff[best*n + pivot]))
                    best = r;
            if (fabs(coeff[best*n + pivot]) < 1e-30) continue;
            if (best != pivot) {
                for (k = 0; k < n; k++) {
                    double t2 = coeff[pivot*n + k];
                    coeff[pivot*n + k] = coeff[best*n + k];
                    coeff[best*n + k] = t2;
                }
                double t2 = rhs[pivot]; rhs[pivot] = rhs[best]; rhs[best] = t2;
            }
            double pv = coeff[pivot*n + pivot];
            for (k = pivot; k < n; k++) coeff[pivot*n + k] /= pv;
            rhs[pivot] /= pv;
            for (r = pivot+1; r < n; r++) {
                double f = coeff[r*n + pivot];
                for (k = pivot; k < n; k++) coeff[r*n + k] -= f * coeff[pivot*n + k];
                rhs[r] -= f * rhs[pivot];
            }
        }
        for (i = n-1; i >= 0; i--) {
            double sv = rhs[i];
            for (k = i+1; k < n; k++) sv -= coeff[i*n + k] * Y[k*n + j];
            Y[i*n + j] = (fabs(coeff[i*n + i]) > 1e-30) ? sv / coeff[i*n + i] : sv;
        }
        free(rhs); free(coeff);
    }
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            double sv = 0.0;
            for (k = 0; k < n; k++) sv += U[i*n + k] * Y[k*n + j];
            temp[i*n + j] = sv;
        }
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            double sv = 0.0;
            for (k = 0; k < n; k++) sv += temp[i*n + k] * U[j*n + k];
            P[i*n + j] = sv;
        }
    for (i = 0; i < n; i++)
        for (j = i+1; j < n; j++) {
            double avg = 0.5 * (P[i*n + j] + P[j*n + i]);
            P[i*n + j] = avg; P[j*n + i] = avg;
        }
    free(T); free(U); free(F); free(temp); free(Y);
    return 0;
}

/* ==============================================================
 * L2 gain (H-infinity norm) via Hamiltonian bisection
 *
 * ||G||_inf = sup_omega sigma_max(G(j*omega))
 *
 * Bounded Real Lemma: ||G||_inf < gamma iff
 * Hamiltonian H = [A, BB^T/gamma^2; -C^T C, -A^T] has no
 * eigenvalues on the imaginary axis.
 *
 * Complexity: O(log(1/tol) * n^3)
 * Reference: Boyd, El Ghaoui, Feron, Balakrishnan (1994) "LMI in
 *   System and Control Theory", Sec 2.7
 * ============================================================== */

double dta_l2_gain(const double* A, const double* B, const double* C,
                    int n, int p, int q) {
    if (!A || !B || !C || n <= 0 || p <= 0 || q <= 0) return INFINITY;
    if (!dta_is_hurwitz(A, n, 1e-10)) return INFINITY;
    int i, j, k;
    double* CTC = safe_alloc((size_t)(n*n) * sizeof(double));
    memset(CTC, 0, (size_t)(n*n) * sizeof(double));
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            for (k = 0; k < q; k++)
                CTC[i*n + j] += C[k*n + i] * C[k*n + j];
    double* Pmat = safe_alloc((size_t)(n*n) * sizeof(double));
    dta_solve_lyapunov(A, n, CTC, Pmat);
    double* BPB = safe_alloc((size_t)(p*p) * sizeof(double));
    double btpb_norm = 0.0;
    for (i = 0; i < p; i++)
        for (j = 0; j < p; j++) {
            double sv = 0.0;
            int k1, k2;
            for (k1 = 0; k1 < n; k1++)
                for (k2 = 0; k2 < n; k2++)
                    sv += B[k1*p + i] * Pmat[k1*n + k2] * B[k2*p + j];
            BPB[i*p + j] = sv;
            if (fabs(sv) > btpb_norm) btpb_norm = fabs(sv);
        }
    double gamma_low = 0.0, gamma_high = 1e6;
    if (btpb_norm > 0) gamma_high = sqrt(btpb_norm) * 2.0;
    int iter;
    for (iter = 0; iter < 30; iter++) {
        double gamma = (gamma_low + gamma_high) / 2.0;
        int m = 2 * n;
        double* Ham = safe_alloc((size_t)(m*m) * sizeof(double));
        memset(Ham, 0, (size_t)(m*m) * sizeof(double));
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                Ham[i*m + j] = A[i*n + j];
                Ham[(n+i)*m + j] = -CTC[i*n + j];
                Ham[(n+i)*m + n + j] = -A[j*n + i];
            }
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                double sv = 0.0;
                for (k = 0; k < p; k++)
                    sv += B[i*p + k] * B[j*p + k];
                Ham[i*m + n + j] = sv / (gamma * gamma);
            }
        double* re = safe_alloc((size_t)m * sizeof(double));
        double* im = safe_alloc((size_t)m * sizeof(double));
        dta_eigenvalues(Ham, m, re, im, 100);
        bool on_imag = false;
        for (i = 0; i < m; i++) {
            if (fabs(re[i]) < 1e-8 && fabs(im[i]) > 1e-10) {
                on_imag = true; break;
            }
        }
        if (on_imag) gamma_low = gamma; else gamma_high = gamma;
        free(Ham); free(re); free(im);
    }
    free(CTC); free(Pmat); free(BPB);
    return (gamma_low + gamma_high) / 2.0;
}
