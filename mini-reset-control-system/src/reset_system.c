/* reset_system.c - Reset Control System Composition and Analysis */
#include "reset_system.h"
#include "reset_math.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

static void hessenberg_reduce(int n, double *A)
{
    for (int k=0; k<n-2; k++) {
        double sigma=0.0;
        for (int i=k+1; i<n; i++) sigma+=A[i*n+k]*A[i*n+k];
        if (sigma<1e-30) continue;
        double x0=A[(k+1)*n+k];
        double alpha=(x0>0.0)?-sqrt(sigma):sqrt(sigma);
        double u1=x0-alpha, beta=1.0/(sigma-x0*alpha);
        for (int j=k; j<n; j++) {
            double dot=u1*A[(k+1)*n+j];
            for (int i=k+2; i<n; i++) dot+=A[i*n+k]*A[i*n+j];
            double s=beta*dot;
            A[(k+1)*n+j]-=s*u1;
            for (int i=k+2; i<n; i++) A[i*n+j]-=s*A[i*n+k];
        }
        for (int i=0; i<n; i++) {
            double dot=u1*A[i*n+(k+1)];
            for (int j=k+2; j<n; j++) dot+=A[j*n+k]*A[i*n+j];
            double s=beta*dot;
            A[i*n+(k+1)]-=s*u1;
            for (int j=k+2; j<n; j++) A[i*n+j]-=s*A[j*n+k];
        }
        A[(k+1)*n+k]=alpha;
        for (int i=k+2; i<n; i++) A[i*n+k]=0.0;
    }
}

static void qr_step_double_shift(int n, double *H, double *er, double *ei, int *p)
{
    int m=n;
    while (m>0) {
        for (int i=m-1; i>=1; i--)
            if (fabs(H[i*n+(i-1)])<1e-12*(fabs(H[(i-1)*n+(i-1)])+fabs(H[i*n+i])))
                H[i*n+(i-1)]=0.0;
        int l=m-1;
        while (l>0&&fabs(H[l*n+(l-1)])>1e-30) l--;
        if (l==m-1) { er[m-1]=H[(m-1)*n+(m-1)]; ei[m-1]=0.0; m--; }
        else if (l==m-2) {
            double a=H[(m-2)*n+(m-2)],b=H[(m-2)*n+(m-1)],c=H[(m-1)*n+(m-2)],d=H[(m-1)*n+(m-1)];
            double tr=a+d, dt=a*d-b*c, disc=tr*tr-4.0*dt;
            if (disc<0.0) {
                er[m-2]=er[m-1]=0.5*tr; ei[m-2]=0.5*sqrt(-disc); ei[m-1]=-ei[m-2];
            } else {
                double sd=sqrt(disc); er[m-2]=0.5*(tr+sd); er[m-1]=0.5*(tr-sd);
                ei[m-2]=ei[m-1]=0.0;
            }
            m-=2;
        } else {
            double s=H[(m-1)*n+(m-1)]+H[(m-2)*n+(m-2)];
            double t=H[(m-1)*n+(m-1)]*H[(m-2)*n+(m-2)]-H[(m-1)*n+(m-2)]*H[(m-2)*n+(m-1)];
            double x=H[l*n+l]*H[l*n+l]+H[l*n+(l+1)]*H[(l+1)*n+l]-s*H[l*n+l]+t;
            double y=H[(l+1)*n+l]*(H[l*n+l]+H[(l+1)*n+(l+1)]-s);
            double z=H[(l+1)*n+l]*H[(l+2)*n+(l+1)];
            for (int k=l; k<m-1; k++) {
                double r=sqrt(x*x+y*y+(k<m-2?z*z:0.0));
                if (r<1e-30) break;
                double v1=x+(x>=0.0?r:-r), v2=y, v3=(k<m-2)?z:0.0, bet=2.0/(v1*v1+v2*v2+v3*v3);
                for (int j=k; j<n; j++) {
                    double d1=bet*(v1*H[k*n+j]+v2*H[(k+1)*n+j]+(k<m-2?v3*H[(k+2)*n+j]:0.0));
                    H[k*n+j]-=d1*v1; H[(k+1)*n+j]-=d1*v2;
                    if(k<m-2) H[(k+2)*n+j]-=d1*v3;
                }
                for (int i=0; i<m&&i<=k+2; i++) {
                    double d2=bet*(v1*H[i*n+k]+v2*H[i*n+(k+1)]+(k<m-2?v3*H[i*n+(k+2)]:0.0));
                    H[i*n+k]-=d2*v1; H[i*n+(k+1)]-=d2*v2;
                    if(k<m-2) H[i*n+(k+2)]-=d2*v3;
                }
                if(k<m-2){x=H[(k+1)*n+k]; y=H[(k+2)*n+k]; z=(k<m-3)?H[(k+3)*n+k]:0.0;}
            }
        }
    }
    *p=m;
}

int reset_eigenvalues(const ResetLinearBase *base, double *er, double *ei)
{
    if (!base||!er||!ei||base->n<=0||!base->A) return 0;
    int n=base->n;
    double *Ac=(double*)malloc((size_t)n*n*sizeof(double));
    if (!Ac) return 0;
    memcpy(Ac,base->A,(size_t)n*n*sizeof(double));
    hessenberg_reduce(n,Ac);
    int max_iter=100*n, iter=0, remaining=n;
    while (iter<max_iter&&remaining>0) {
        int old_r=remaining;
        qr_step_double_shift(n,Ac,er,ei,&remaining);
        if (remaining>=old_r) break;
        iter++;
    }
    if (remaining>0) for (int i=0; i<remaining; i++) { er[i]=0.0; ei[i]=0.0; }
    free(Ac);
    return n;
}

bool reset_is_hurwitz(const ResetLinearBase *base)
{
    if (!base||!base->A||base->n<=0) return false;
    int n=base->n;
    double *er=(double*)malloc((size_t)n*sizeof(double));
    double *ei=(double*)malloc((size_t)n*sizeof(double));
    if (!er||!ei) { free(er); free(ei); return false; }
    int nf=reset_eigenvalues(base,er,ei);
    bool h=true;
    for (int i=0; i<nf; i++) if (er[i]>=-1e-10) { h=false; break; }
    free(er); free(ei);
    return h;
}

double reset_h2_norm(const ResetLinearBase *base)
{
    if (!base||!base->A||base->n<=0||!reset_is_hurwitz(base)) return INFINITY;
    int n=base->n, m=base->m, p=base->p;
    double *Q=(double*)calloc((size_t)n*n,sizeof(double));
    if (!Q) return INFINITY;
    for (int i=0; i<n; i++) for (int j=0; j<n; j++)
        for (int k=0; k<m; k++) Q[i*n+j]+=base->B[i*m+k]*base->B[j*m+k];
    double *P=(double*)calloc((size_t)n*n,sizeof(double));
    if (!P) { free(Q); return INFINITY; }
    if (mat_lyapunov_solve(n,base->A,Q,P)!=0) { free(Q); free(P); return INFINITY; }
    free(Q);
    double h2sq=0.0;
    for (int i=0; i<p; i++) for (int k=0; k<n; k++) {
        double cp=0.0;
        for (int l=0; l<n; l++) cp+=base->C[i*n+l]*P[l*n+k];
        h2sq+=cp*base->C[i*n+k];
    }
    if (base->D) for (int i=0; i<p; i++) for (int j=0; j<m; j++)
        h2sq+=base->D[i*m+j]*base->D[i*m+j];
    free(P);
    return sqrt(h2sq);
}

double reset_hinf_norm(const ResetLinearBase *base, int mi, double tol)
{
    if (!base||!base->A||base->n<=0) return INFINITY;
    (void)mi; (void)tol;
    double peak=0.0;
    for (int k=0; k<600; k++) {
        double w=1e-4*pow(1e4/1e-4,(double)k/599.0);
        double sv=mat_sigma_max(base->n,base->m,base->p,base->A,base->B,base->C,base->D,w,30,1e-6);
        if (sv>peak) peak=sv;
    }
    return peak;
}

double* reset_closed_loop_A(const ResetLinearBase *plant,
                             const ResetLinearBase *controller, int *n_cl)
{
    if (!plant||!controller||!n_cl) return NULL;
    int np=plant->n, nc=controller->n, nt=np+nc;
    *n_cl=nt;
    double *Ac=(double*)calloc((size_t)nt*nt,sizeof(double));
    if (!Ac) return NULL;
    for (int i=0; i<np; i++) for (int j=0; j<np; j++) Ac[i*nt+j]=plant->A[i*np+j];
    if (plant->B&&controller->C&&plant->m==controller->p)
        for (int i=0; i<np; i++) for (int j=0; j<nc; j++)
            for (int k=0; k<plant->m; k++)
                Ac[i*nt+(np+j)]+=plant->B[i*plant->m+k]*controller->C[k*nc+j];
    if (controller->B&&plant->C&&controller->m==plant->p)
        for (int i=0; i<nc; i++) for (int j=0; j<np; j++)
            for (int k=0; k<controller->m; k++)
                Ac[(np+i)*nt+j]-=controller->B[i*controller->m+k]*plant->C[k*np+j];
    for (int i=0; i<nc; i++) for (int j=0; j<nc; j++) Ac[(np+i)*nt+(np+j)]=controller->A[i*nc+j];
    return Ac;
}

bool reset_check_hbeta_stability(const ResetLinearBase *plant,
                                  const ResetSystem *controller)
{
    if (!plant||!controller||!controller->flow) return false;
    int ncl; double *Acl=reset_closed_loop_A(plant,controller->flow,&ncl);
    if (!Acl) return false;
    ResetLinearBase b; memset(&b,0,sizeof(b)); b.n=ncl; b.A=Acl;
    bool stable=reset_is_hurwitz(&b);
    free(Acl);
    return stable;
}

ResetFeedbackLoop* reset_feedback_create(const ResetSystem *c, const ResetLinearBase *p)
{
    if (!c||!p) return NULL;
    ResetFeedbackLoop *l=(ResetFeedbackLoop*)calloc(1,sizeof(ResetFeedbackLoop));
    if (!l) return NULL;
    l->plant=reset_base_clone(p); l->controller=(ResetSystem*)c;
    l->np=p->n; l->nc=c->nc;
    l->xp=(double*)calloc((size_t)p->n,sizeof(double));
    l->xc=c->x_c; l->dt=0.001; l->is_initialized=true;
    if (!l->plant||!l->xp) { reset_feedback_free(l); return NULL; }
    return l;
}

void reset_feedback_free(ResetFeedbackLoop *l)
{
    if (!l) return;
    if (l->plant) reset_base_free(l->plant);
    free(l->xp); free(l);
}

double reset_feedback_step(ResetFeedbackLoop *l, double dt, double r)
{
    if (!l||!l->plant||!l->controller) return 0.0;
    l->r=r; l->e=r-l->y;
    reset_check_and_apply(l->controller,l->e,l->controller->e_prev);
    double u=0.0; ResetLinearBase *cf=l->controller->flow;
    if (cf&&cf->C&&l->xc) for(int j=0;j<l->nc;j++) u+=cf->C[j]*l->xc[j];
    if (cf&&cf->D) u+=cf->D[0]*l->e;
    l->u=u;
    if (cf&&cf->A&&cf->B&&l->xc) for(int i=0;i<l->nc;i++) {
        double dx=0.0;
        for(int j=0;j<l->nc;j++) dx+=cf->A[i*l->nc+j]*l->xc[j];
        dx+=cf->B[i]*l->e; l->xc[i]+=dx*dt;
    }
    ResetLinearBase *pf=l->plant;
    if (pf&&pf->A&&pf->B&&l->xp) for(int i=0;i<l->np;i++) {
        double dx=0.0;
        for(int j=0;j<l->np;j++) dx+=pf->A[i*l->np+j]*l->xp[j];
        dx+=pf->B[i]*l->u; l->xp[i]+=dx*dt;
    }
    l->y=0.0;
    if (pf&&pf->C&&l->xp) for(int j=0;j<l->np;j++) l->y+=pf->C[j]*l->xp[j];
    if (pf&&pf->D) l->y+=pf->D[0]*l->u;
    l->controller->t+=dt; l->controller->e_prev=l->e;
    return l->y;
}

int reset_feedback_simulate(ResetFeedbackLoop *l, const double *ts,
                             const double *rv, double *yo, int ns)
{
    if (!l||!ts||!rv||!yo||ns<=0) return 0;
    int cnt=0;
    for(int k=0;k<ns;k++) {
        double dt=(k==0)?ts[0]:ts[k]-ts[k-1]; if(dt<=0.0) dt=0.001;
        yo[k]=reset_feedback_step(l,dt,rv[k]);
        const ResetInterval *ri=reset_get_interval_stats(l->controller);
        if(ri) cnt=ri->n_resets;
    }
    return cnt;
}

double reset_closed_loop_freqresp(const ResetLinearBase *p, const ResetSystem *c, double w)
{
    if (!p||!c||!c->flow) return 0.0;
    if (p->p==1&&p->m==1&&c->flow->p==1&&c->flow->m==1) {
        double Gpr,Gpi,Gcr,Gci;
        mat_transfer_func(p->n,1,1,p->A,p->B,p->C,p->D,0.0,w,&Gpr,&Gpi);
        mat_transfer_func(c->flow->n,1,1,c->flow->A,c->flow->B,c->flow->C,c->flow->D,0.0,w,&Gcr,&Gci);
        double Lr=Gpr*Gcr-Gpi*Gci, Li=Gpr*Gci+Gpi*Gcr;
        double dr=1.0+Lr, di=Li, den=dr*dr+di*di;
        if(den<1e-30) return 1e6;
        double Tr=(Lr*dr+Li*di)/den, Ti=(Li*dr-Lr*di)/den;
        return sqrt(Tr*Tr+Ti*Ti);
    }
    return 0.0;
}

double reset_sensitivity(const ResetLinearBase *p, const ResetSystem *c, double w)
{
    if (!p||!c||!c->flow) return 1.0;
    if (p->p==1&&p->m==1) {
        double Gpr,Gpi,Gcr,Gci;
        mat_transfer_func(p->n,1,1,p->A,p->B,p->C,p->D,0.0,w,&Gpr,&Gpi);
        mat_transfer_func(c->flow->n,1,1,c->flow->A,c->flow->B,c->flow->C,c->flow->D,0.0,w,&Gcr,&Gci);
        double Lr=Gpr*Gcr-Gpi*Gci, Li=Gpr*Gci+Gpi*Gcr;
        return 1.0/sqrt((1.0+Lr)*(1.0+Lr)+Li*Li);
    }
    return 1.0;
}

double reset_complementary_sensitivity(const ResetLinearBase *p, const ResetSystem *c, double w)
{ return reset_closed_loop_freqresp(p,c,w); }

ResetSystem* reset_series_connection(const ResetLinearBase *s1, const ResetSystem *s2)
{
    if (!s1||!s2) return NULL;
    ResetSystem *r=reset_sys_create(s2->nc);
    if (!r) return NULL;
    reset_base_free(r->flow); r->flow=reset_base_clone(s2->flow);
    if (s2->jump&&r->jump) {
        int n2=s2->jump->nc; free(r->jump->Ar); free(r->jump->Br);
        r->jump->nc=n2;
        r->jump->Ar=(double*)malloc((size_t)n2*n2*sizeof(double));
        r->jump->Br=(double*)malloc((size_t)n2*sizeof(double));
        if(r->jump->Ar&&s2->jump->Ar) memcpy(r->jump->Ar,s2->jump->Ar,(size_t)n2*n2*sizeof(double));
    }
    return r;
}

ResetLinearBase* reset_feedback_connection(const ResetLinearBase *p, const ResetSystem *c)
{
    if (!p||!c||!c->flow) return NULL;
    int ncl; double *Acl=reset_closed_loop_A(p,c->flow,&ncl);
    if (!Acl) return NULL;
    ResetLinearBase *fb=reset_base_create(ncl,p->m,p->p);
    if (!fb) { free(Acl); return NULL; }
    memcpy(fb->A,Acl,(size_t)ncl*ncl*sizeof(double));
    free(Acl);
    return fb;
}
