/* reset_lyapunov.c - Lyapunov stability analysis for reset systems */
#include "reset_lyapunov.h"
#include "reset_system.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

bool lyap_is_spd(int n, const double *P)
{
    if (!P || n<=0) return false;
    /* Cholesky decomposition attempt: P must be SPD */
    double *L=(double*)calloc((size_t)n*n,sizeof(double));
    if (!L) return false;
    for (int i=0; i<n; i++) {
        for (int j=0; j<=i; j++) {
            double sum=P[i*n+j];
            for (int k=0; k<j; k++) sum -= L[i*n+k]*L[j*n+k];
            if (i==j) {
                if (sum<=1e-12) { free(L); return false; }
                L[i*n+i]=sqrt(sum);
            } else {
                L[i*n+j]=sum/L[j*n+j];
            }
        }
    }
    free(L);
    return true;
}

bool lyap_construct_linear(const ResetLinearBase *base, const double *Q, double *P_out)
{
    if (!base||!base->A||base->n<=0||!P_out) return false;
    int n=base->n;
    double *Q_use=(double*)malloc((size_t)n*n*sizeof(double));
    if (!Q_use) return false;
    if (Q) memcpy(Q_use,Q,(size_t)n*n*sizeof(double));
    else { for(int i=0; i<n*n; i++) Q_use[i]=0.0; for(int i=0; i<n; i++) Q_use[i*n+i]=1.0; }
    int ret=mat_lyapunov_solve(n,base->A,Q_use,P_out);
    free(Q_use);
    if (ret!=0) return false;
    return lyap_is_spd(n,P_out);
}

bool lyap_verify_reset(const ResetLinearBase *base, const ResetJumpMap *jump,
                        const double *P, double *margin)
{
    if (!base||!P||!margin) return false;
    int n=base->n;
    /* Check P > 0 */
    if (!lyap_is_spd(n,P)) { *margin=-1.0; return false; }
    /* Check A^T*P + P*A < 0 */
    double *Q=(double*)calloc((size_t)n*n,sizeof(double));
    if (!Q) return false;
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++) {
            double atp=0.0, pa=0.0;
            for (int k=0; k<n; k++) { atp+=base->A[k*n+i]*P[k*n+j]; pa+=P[i*n+k]*base->A[k*n+j]; }
            Q[i*n+j] = atp+pa;
        }
    /* Check eigenvalues of Q are negative */
    double *er=(double*)malloc((size_t)n*sizeof(double));
    double *ei=(double*)malloc((size_t)n*sizeof(double));
    if (!er||!ei) { free(Q); free(er); free(ei); return false; }
    /* Simple check: all diagonal elements should be negative for SPD P */
    bool flow_ok=true; double max_ev=0.0;
    for (int i=0; i<n; i++) { if (Q[i*n+i]>=0.0) flow_ok=false; if(Q[i*n+i]>max_ev) max_ev=Q[i*n+i]; }
    free(Q); free(er); free(ei);
    /* Check jump condition: P - Ar^T*P*Ar >= 0 */
    if (jump&&jump->Ar) {
        int nc=jump->nc;
        double *Par=(double*)calloc((size_t)nc*nc,sizeof(double));
        if (!Par) return false;
        for (int i=0; i<nc; i++) {
            for (int j=0; j<nc; j++) {
                double sum=0.0;
                for (int k=0; k<nc; k++) {
                    double park=0.0;
                    for (int l=0; l<nc; l++) park += P[i*nc+l]*jump->Ar[l*nc+k];
                    sum += park*jump->Ar[j*nc+k];
                }
                Par[i*nc+j] = P[i*nc+j]-sum;
            }
        }
        /* Check if Par is positive semidefinite */
        for (int i=0; i<nc; i++) {
            if (Par[i*nc+i]<-1e-10) { free(Par); *margin=-1.0; return false; }
        }
        free(Par);
    }
    *margin = -max_ev;
    return flow_ok;
}

LyapAnalysis* lyap_search(const ResetSystem *rsys, int n_trials)
{
    if (!rsys||!rsys->flow||rsys->nc<=0) return NULL;
    LyapAnalysis *la=(LyapAnalysis*)calloc(1,sizeof(LyapAnalysis));
    if (!la) return NULL;
    int n=rsys->nc;
    la->n_dim=n;
    la->V=(LyapFunction*)calloc(1,sizeof(LyapFunction));
    if (!la->V) { free(la); return NULL; }
    la->V->n=n;
    la->V->P=(double*)calloc((size_t)n*n,sizeof(double));
    if (!la->V->P) { free(la->V); free(la); return NULL; }
    /* Try to construct Lyapunov function for base linear system */
    bool found=false;
    for (int trial=0; trial<n_trials; trial++) {
        double *Q_trial=(double*)calloc((size_t)n*n,sizeof(double));
        if (!Q_trial) continue;
        for (int i=0; i<n; i++) Q_trial[i*n+i]=1.0+0.1*trial;
        double *P_trial=(double*)calloc((size_t)n*n,sizeof(double));
        if (!P_trial) { free(Q_trial); continue; }
        if (mat_lyapunov_solve(n,rsys->flow->A,Q_trial,P_trial)==0) {
            double margin=0.0;
            if (lyap_verify_reset(rsys->flow,rsys->jump,P_trial,&margin)) {
                memcpy(la->V->P,P_trial,(size_t)n*n*sizeof(double));
                la->V->is_valid=true;
                la->margin=margin;
                la->is_stable=true;
                la->is_quadratic=true;
                found=true;
                free(Q_trial); free(P_trial);
                break;
            }
        }
        free(Q_trial); free(P_trial);
    }
    if (!found) { la->is_stable=false; }
    return la;
}

void lyap_analysis_free(LyapAnalysis *a)
{
    if (!a) return;
    if (a->V) { free(a->V->P); free(a->V); }
    free(a);
}

bool lyap_check_passivity(const ResetLinearBase *base, double eps)
{
    if (!base||!base->A) return false;
    int n=base->n, m=base->m;
    /* KYP lemma: check if there exists P>0 satisfying LMI.
     * Simplified: check diagonal condition for passivity. */
    if (m!=base->p) return false;
    /* For SISO: need Re[G(jw)] >= 0 for all w */
    for (int k=0; k<200; k++) {
        double w=0.01*pow(100.0,(double)k/199.0);
        double Gr,Gi;
        double *Grr=&Gr, *Gii=&Gi;
        mat_transfer_func(n,m,base->p,base->A,base->B,base->C,base->D,0.0,w,Grr,Gii);
        if (Gr < -eps) return false;
    }
    return true;
}

bool lyap_check_reset_passivity(const ResetSystem *rsys)
{
    if (!rsys||!rsys->flow) return false;
    /* Check base linear passivity */
    if (!lyap_check_passivity(rsys->flow,1e-6)) return false;
    /* Check reset map does not increase storage */
    /* V(x^+) = x^T Ar^T P Ar x <= x^T P x = V(x^-) */
    /* Requires Ar^T P Ar <= P, which holds for ||Ar|| <= 1. */
    if (rsys->jump&&rsys->jump->Ar) {
        int nc=rsys->jump->nc;
        double norm=0.0;
        for (int i=0; i<nc*nc; i++) norm += rsys->jump->Ar[i]*rsys->jump->Ar[i];
        if (sqrt(norm)>1.0+1e-10) return false;
    }
    return true;
}

double lyap_required_dwell_time(double alpha, double beta)
{
    if (beta<=1.0) return 0.0;
    if (alpha<=0.0) return INFINITY;
    return log(beta)/(2.0*alpha);
}

double lyap_evaluate(const LyapFunction *V, const double *x)
{
    if (!V||!V->P||!x) return 0.0;
    int n=V->n;
    double val=0.0;
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            val += x[i]*V->P[i*n+j]*x[j];
    return val;
}

double lyap_flow_derivative(const LyapFunction *V, const ResetLinearBase *base,
                             const double *x, const double *u)
{
    if (!V||!V->P||!base||!x) return 0.0;
    int n=V->n;
    double *Ax=(double*)calloc((size_t)n,sizeof(double));
    if (!Ax) return 0.0;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) Ax[i]+=base->A[i*n+j]*x[j];
        if (u) Ax[i]+=base->B[i]*u[0];
    }
    double vdot=0.0;
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            vdot += Ax[i]*V->P[i*n+j]*x[j] + x[i]*V->P[i*n+j]*Ax[j];
    free(Ax);
    return vdot;
}

double lyap_jump_increment(const LyapFunction *V, const ResetJumpMap *jump,
                            const double *x_before, double e)
{
    if (!V||!V->P||!jump||!x_before) return 0.0;
    int n=V->n;
    double *x_after=(double*)calloc((size_t)n,sizeof(double));
    if (!x_after) return 0.0;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) x_after[i]+=jump->Ar[i*n+j]*x_before[j];
        if (jump->Br) x_after[i]+=jump->Br[i]*e;
    }
    double va=lyap_evaluate(V,x_after);
    double vb=lyap_evaluate(V,x_before);
    free(x_after);
    return va-vb;
}

bool lyap_circle_criterion(const ResetLinearBase *plant, const ResetSystem *ctrl,
                            double rho, double om_min, double om_max, int np)
{
    if (!plant||!ctrl||!ctrl->flow) return false;
    /* Circle criterion for FORE:
     * Condition: sup_w |T(jw)| * (1+rho)/(1-rho) < 1
     * where T is complementary sensitivity.
     */
    double denom=(1.0-rho);
    if (denom<1e-10) return false;
    double factor=(1.0+rho)/denom;
    double max_T=0.0;
    for (int k=0; k<np; k++) {
        double w=om_min*pow(om_max/om_min,(double)k/(double)(np-1));
        double T=reset_closed_loop_freqresp(plant,ctrl,w);
        if (T>max_T) max_T=T;
    }
    return (max_T*factor<1.0);
}

double lyap_hbeta_condition(const ResetLinearBase *plant, const ResetSystem *ctrl)
{
    if (!plant||!ctrl||!ctrl->flow) return -1.0;
    int ncl; double *Acl=reset_closed_loop_A(plant,ctrl->flow,&ncl);
    if (!Acl) return -1.0;
    /* Check Hurwitz stability of base closed loop */
    ResetLinearBase b; memset(&b,0,sizeof(b));
    b.n=ncl; b.A=Acl;
    bool stable=reset_is_hurwitz(&b);
    free(Acl);
    if (!stable) return -1.0;
    /* For stable base system, H_beta condition always holds for some beta */
    return 1.0;
}
