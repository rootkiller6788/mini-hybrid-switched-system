/*
 * impulsive_flow.c -- Continuous flow dynamics implementations
 *
 * Provides vector field evaluations for canonical dynamical systems
 * used as testbeds in impulsive control and synchronization.
 *
 * Each system implements a distinct ImpVectorField callback.
 * Systems covered: linear, Lorenz, Chua, Duffing, Rossler,
 * Van der Pol, FitzHugh-Nagumo, pendulum, harmonic oscillator.
 *
 * References:
 *   Lorenz (1963) J. Atmos. Sci. 20:130-141
 *   Chua (1992) "The Genesis of Chua's Circuit"
 *   FitzHugh (1961) Biophys. J. 1:445-466
 *   Strogatz (2015) "Nonlinear Dynamics and Chaos"
 */
#include "impulsive_flow.h"
#include <stdlib.h>
#include <math.h>

/* ---- Linear Flow ---- */

ImpFlowLinear* affine_flow_create(const double *A, int n)
{
    if (!A || n < 1 || n > IMP_MAX_DIM) return NULL;
    ImpFlowLinear *flow = (ImpFlowLinear*)calloc(1, sizeof(ImpFlowLinear));
    if (!flow) return NULL;
    flow->A = (double*)malloc((size_t)n * n * sizeof(double));
    if (!flow->A) { free(flow); return NULL; }
    for (int i = 0; i < n * n; i++) flow->A[i] = A[i];
    flow->n = n;
    return flow;
}

void affine_flow_free(ImpFlowLinear *flow)
{
    if (!flow) return;
    free(flow->A); free(flow);
}

int affine_flow_eval(double t, const double *x, int n, double *dxdt, void *ctx)
{
    (void)t;
    ImpFlowLinear *flow = (ImpFlowLinear*)ctx;
    if (!flow || !x || !dxdt || n != flow->n) return -1;
    for (int i = 0; i < n; i++) {
        dxdt[i] = 0.0;
        for (int j = 0; j < n; j++)
            dxdt[i] += flow->A[i * n + j] * x[j];
    }
    return 0;
}

/* ---- Lorenz System ---- */

ImpFlowLorenz* imp_flow_lorenz_create(double sigma, double rho, double beta)
{
    if (!isfinite(sigma) || !isfinite(rho) || !isfinite(beta)) return NULL;
    if (sigma <= 0.0 || beta <= 0.0) return NULL;
    ImpFlowLorenz *lrz = (ImpFlowLorenz*)calloc(1, sizeof(ImpFlowLorenz));
    if (!lrz) return NULL;
    lrz->sigma = sigma; lrz->rho = rho; lrz->beta = beta;
    return lrz;
}

void imp_flow_lorenz_free(ImpFlowLorenz *lrz) { free(lrz); }

int imp_flow_lorenz_eval(double t, const double *x, int n,
                          double *dxdt, void *ctx)
{
    (void)t;
    ImpFlowLorenz *lrz = (ImpFlowLorenz*)ctx;
    if (!lrz || !x || !dxdt || n != 3) return -1;
    double s = lrz->sigma, r = lrz->rho, b = lrz->beta;
    dxdt[0] = s * (x[1] - x[0]);
    dxdt[1] = x[0] * (r - x[2]) - x[1];
    dxdt[2] = x[0] * x[1] - b * x[2];
    return 0;
}

/* ---- Chua's Circuit ---- */

ImpFlowChua* imp_flow_chua_create(double alpha, double beta,
                                    double m0, double m1)
{
    if (!isfinite(alpha) || !isfinite(beta) ||
        !isfinite(m0) || !isfinite(m1)) return NULL;
    if (alpha <= 0.0 || beta <= 0.0) return NULL;
    ImpFlowChua *chua = (ImpFlowChua*)calloc(1, sizeof(ImpFlowChua));
    if (!chua) return NULL;
    chua->alpha = alpha; chua->beta = beta;
    chua->m0 = m0; chua->m1 = m1;
    return chua;
}

void imp_flow_chua_free(ImpFlowChua *chua) { free(chua); }

static double chua_phi(double x, double m0, double m1)
{
    return m1 * x + 0.5 * (m0 - m1) * (fabs(x + 1.0) - fabs(x - 1.0));
}

int imp_flow_chua_eval(double t, const double *x, int n,
                        double *dxdt, void *ctx)
{
    (void)t;
    ImpFlowChua *ch = (ImpFlowChua*)ctx;
    if (!ch || !x || !dxdt || n != 3) return -1;
    double a = ch->alpha, b = ch->beta;
    double phi_x0 = chua_phi(x[0], ch->m0, ch->m1);
    dxdt[0] = a * (x[1] - x[0] - phi_x0);
    dxdt[1] = x[0] - x[1] + x[2];
    dxdt[2] = -b * x[1];
    return 0;
}

/* ---- Duffing Oscillator ---- */

ImpFlowDuffing* imp_flow_duffing_create(double delta, double alpha,
    double beta, double gamma, double omega)
{
    if (!isfinite(delta) || !isfinite(alpha) || !isfinite(beta) ||
        !isfinite(gamma) || !isfinite(omega)) return NULL;
    ImpFlowDuffing *duf = (ImpFlowDuffing*)calloc(1, sizeof(ImpFlowDuffing));
    if (!duf) return NULL;
    duf->delta = delta; duf->alpha = alpha;
    duf->beta = beta; duf->gamma = gamma; duf->omega = omega;
    return duf;
}

void imp_flow_duffing_free(ImpFlowDuffing *duf) { free(duf); }

int imp_flow_duffing_eval(double t, const double *x, int n,
                           double *dxdt, void *ctx)
{
    ImpFlowDuffing *duf = (ImpFlowDuffing*)ctx;
    if (!duf || !x || !dxdt || n != 2) return -1;
    double d = duf->delta, a = duf->alpha, b = duf->beta;
    double g = duf->gamma, w = duf->omega;
    dxdt[0] = x[1];
    dxdt[1] = -d * x[1] + a * x[0] - b * x[0] * x[0] * x[0]
              + g * cos(w * t);
    return 0;
}

/* ---- Rossler Attractor ---- */

ImpFlowRossler* imp_flow_rossler_create(double a, double b, double c)
{
    if (!isfinite(a) || !isfinite(b) || !isfinite(c)) return NULL;
    ImpFlowRossler *r = (ImpFlowRossler*)calloc(1, sizeof(ImpFlowRossler));
    if (!r) return NULL;
    r->a = a; r->b = b; r->c = c;
    return r;
}

void imp_flow_rossler_free(ImpFlowRossler *r) { free(r); }

int imp_flow_rossler_eval(double t, const double *x, int n,
                           double *dxdt, void *ctx)
{
    (void)t;
    ImpFlowRossler *rr = (ImpFlowRossler*)ctx;
    if (!rr || !x || !dxdt || n != 3) return -1;
    dxdt[0] = -x[1] - x[2];
    dxdt[1] = x[0] + rr->a * x[1];
    dxdt[2] = rr->b + x[2] * (x[0] - rr->c);
    return 0;
}

/* ---- Van der Pol Oscillator ---- */

ImpFlowVanDerPol* imp_flow_vdp_create(double mu)
{
    if (!isfinite(mu)) return NULL;
    ImpFlowVanDerPol *vdp = (ImpFlowVanDerPol*)calloc(1, sizeof(ImpFlowVanDerPol));
    if (!vdp) return NULL;
    vdp->mu = mu;
    return vdp;
}

void imp_flow_vdp_free(ImpFlowVanDerPol *vdp) { free(vdp); }

int imp_flow_vdp_eval(double t, const double *x, int n,
                       double *dxdt, void *ctx)
{
    (void)t;
    ImpFlowVanDerPol *vd = (ImpFlowVanDerPol*)ctx;
    if (!vd || !x || !dxdt || n != 2) return -1;
    dxdt[0] = x[1];
    dxdt[1] = vd->mu * (1.0 - x[0] * x[0]) * x[1] - x[0];
    return 0;
}

/* ---- FitzHugh-Nagumo Neuron Model ---- */

ImpFlowFHN* imp_flow_fhn_create(double epsilon, double a, double b,
                                 double I_ext)
{
    if (!isfinite(epsilon) || !isfinite(a) ||
        !isfinite(b) || !isfinite(I_ext)) return NULL;
    ImpFlowFHN *fhn = (ImpFlowFHN*)calloc(1, sizeof(ImpFlowFHN));
    if (!fhn) return NULL;
    fhn->epsilon = epsilon; fhn->a = a;
    fhn->b = b; fhn->I_ext = I_ext;
    return fhn;
}

void imp_flow_fhn_free(ImpFlowFHN *fhn) { free(fhn); }

int imp_flow_fhn_eval(double t, const double *x, int n,
                       double *dxdt, void *ctx)
{
    (void)t;
    ImpFlowFHN *fn = (ImpFlowFHN*)ctx;
    if (!fn || !x || !dxdt || n != 2) return -1;
    double eps = fn->epsilon, a = fn->a, b = fn->b, I = fn->I_ext;
    dxdt[0] = x[0] - x[0] * x[0] * x[0] / 3.0 - x[1] + I;
    dxdt[1] = eps * (x[0] + a - b * x[1]);
    return 0;
}

/* ---- Damped Pendulum ---- */

ImpFlowPendulum* imp_flow_pendulum_create(double g_over_L, double b)
{
    if (!isfinite(g_over_L) || !isfinite(b)) return NULL;
    ImpFlowPendulum *pen = (ImpFlowPendulum*)calloc(1, sizeof(ImpFlowPendulum));
    if (!pen) return NULL;
    pen->g_over_L = g_over_L; pen->b = b;
    return pen;
}

void imp_flow_pendulum_free(ImpFlowPendulum *pen) { free(pen); }

int imp_flow_pendulum_eval(double t, const double *x, int n,
                            double *dxdt, void *ctx)
{
    (void)t;
    ImpFlowPendulum *pn = (ImpFlowPendulum*)ctx;
    if (!pn || !x || !dxdt || n != 2) return -1;
    dxdt[0] = x[1];
    dxdt[1] = -pn->g_over_L * sin(x[0]) - pn->b * x[1];
    return 0;
}

/* ---- Harmonic Oscillator ---- */

ImpFlowHarmonic* imp_flow_harmonic_create(double omega0)
{
    if (!isfinite(omega0) || omega0 <= 0.0) return NULL;
    ImpFlowHarmonic *h = (ImpFlowHarmonic*)calloc(1, sizeof(ImpFlowHarmonic));
    if (!h) return NULL;
    h->omega0_sq = omega0 * omega0;
    return h;
}

void imp_flow_harmonic_free(ImpFlowHarmonic *h) { free(h); }

int imp_flow_harmonic_eval(double t, const double *x, int n,
                            double *dxdt, void *ctx)
{
    (void)t;
    ImpFlowHarmonic *ho = (ImpFlowHarmonic*)ctx;
    if (!ho || !x || !dxdt || n != 2) return -1;
    dxdt[0] = x[1];
    dxdt[1] = -ho->omega0_sq * x[0];
    return 0;
}

/* ---- Generic Flow Wrapper ---- */

ImpFlowGeneric* imp_flow_generic_create(ImpVectorField vf, void *params, int n)
{
    if (!vf || n < 1 || n > IMP_MAX_DIM) return NULL;
    ImpFlowGeneric *gen = (ImpFlowGeneric*)calloc(1, sizeof(ImpFlowGeneric));
    if (!gen) return NULL;
    gen->vf = vf; gen->params = params; gen->n = n;
    return gen;
}

void imp_flow_generic_free(ImpFlowGeneric *gen) { free(gen); }

int imp_flow_generic_eval(double t, const double *x, int n,
                           double *dxdt, void *ctx)
{
    (void)ctx;
    ImpFlowGeneric *g = (ImpFlowGeneric*)ctx;
    if (!g || !x || !dxdt || n != g->n) return -1;
    return g->vf(t, x, n, dxdt, g->params);
}
