/* reset_math.c - Matrix operations and equation solvers */
#include "reset_math.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

void mat_vec_mul(int m, int n, const double *A, const double *x, double *y, double alpha, double beta)
{
    if (!A || !x || !y || m<=0 || n<=0) return;
    for (int i=0; i<m; i++) y[i] *= beta;
    for (int i=0; i<m; i++) {
        double s=0.0;
        for (int j=0; j<n; j++) s += A[i*n+j]*x[j];
        y[i] += alpha*s;
    }
}

void mat_mat_mul(int m, int k, int n, const double *A, const double *B, double *C, double alpha, double beta)
{
    if (!A||!B||!C||m<=0||k<=0||n<=0) return;
    for (int i=0; i<m; i++) for (int j=0; j<n; j++) C[i*n+j] *= beta;
    for (int i=0; i<m; i++) for (int j=0; j<n; j++) {
        double s=0.0;
        for (int l=0; l<k; l++) s += A[i*k+l]*B[l*n+j];
        C[i*n+j] += alpha*s;
    }
}

void mat_transpose(int m, int n, const double *A, double *B)
{
    if (!A||!B||m<=0||n<=0) return;
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            B[j*m+i] = A[i*n+j];
}

int mat_invert(int n, const double *A, double *A_inv)
{
    if (!A||!A_inv||n<=0) return -1;
    double *LU=(double*)malloc((size_t)n*n*sizeof(double));
    if (!LU) return -1;
    memcpy(LU,A,(size_t)n*n*sizeof(double));
    int *piv=(int*)malloc((size_t)n*sizeof(int));
    if (!piv) { free(LU); return -1; }
    for (int i=0; i<n; i++) piv[i]=i;
    for (int col=0; col<n; col++) {
        double mv=0.0; int mr=col;
        for (int r=col; r<n; r++)
            if (fabs(LU[r*n+col])>mv) { mv=fabs(LU[r*n+col]); mr=r; }
        if (mv<1e-15) { free(LU); free(piv); return -1; }
        if (mr!=col) {
            for (int j=0; j<n; j++) {
                double t=LU[col*n+j]; LU[col*n+j]=LU[mr*n+j]; LU[mr*n+j]=t;
            }
            int t=piv[col]; piv[col]=piv[mr]; piv[mr]=t;
        }
        for (int r=col+1; r<n; r++) {
            double f=LU[r*n+col]/LU[col*n+col]; LU[r*n+col]=f;
            for (int j=col+1; j<n; j++) LU[r*n+j] -= f*LU[col*n+j];
        }
    }
    for (int k=0; k<n; k++) {
        double *b=(double*)calloc((size_t)n,sizeof(double));
        if (!b) { free(LU); free(piv); return -1; }
        b[piv[k]]=1.0;
        for (int i=1; i<n; i++) {
            double s=0.0;
            for (int j=0; j<i; j++) s += LU[i*n+j]*b[j];
            b[i] -= s;
        }
        for (int i=n-1; i>=0; i--) {
            double s=0.0;
            for (int j=i+1; j<n; j++) s += LU[i*n+j]*b[j];
            b[i] = (b[i]-s)/LU[i*n+i];
        }
        for (int i=0; i<n; i++) A_inv[i*n+k]=b[i];
        free(b);
    }
    free(LU); free(piv);
    return 0;
}

int mat_solve(int n, const double *A, double *bx)
{
    if (!A||!bx||n<=0) return -1;
    double *LU=(double*)malloc((size_t)n*n*sizeof(double));
    if (!LU) return -1;
    memcpy(LU,A,(size_t)n*n*sizeof(double));
    for (int col=0; col<n; col++) {
        double mv=0.0; int mr=col;
        for (int r=col; r<n; r++)
            if (fabs(LU[r*n+col])>mv) { mv=fabs(LU[r*n+col]); mr=r; }
        if (mv<1e-15) { free(LU); return -1; }
        if (mr!=col) {
            for (int j=0; j<n; j++) {
                double t=LU[col*n+j]; LU[col*n+j]=LU[mr*n+j]; LU[mr*n+j]=t;
            }
            double t=bx[col]; bx[col]=bx[mr]; bx[mr]=t;
        }
        for (int r=col+1; r<n; r++) {
            double f=LU[r*n+col]/LU[col*n+col];
            for (int j=col+1; j<n; j++) LU[r*n+j] -= f*LU[col*n+j];
            bx[r] -= f*bx[col];
        }
    }
    for (int i=n-1; i>=0; i--) {
        double s=0.0;
        for (int j=i+1; j<n; j++) s += LU[i*n+j]*bx[j];
        bx[i] = (bx[i]-s)/LU[i*n+i];
    }
    free(LU);
    return 0;
}

double mat_determinant(int n, const double *A)
{
    if (!A||n<=0) return 0.0;
    double *LU=(double*)malloc((size_t)n*n*sizeof(double));
    if (!LU) return 0.0;
    memcpy(LU,A,(size_t)n*n*sizeof(double));
    double det=1.0; int sign=1;
    for (int col=0; col<n; col++) {
        double mv=0.0; int mr=col;
        for (int r=col; r<n; r++) if (fabs(LU[r*n+col])>mv) { mv=fabs(LU[r*n+col]); mr=r; }
        if (mv<1e-15) { free(LU); return 0.0; }
        if (mr!=col) {
            for (int j=0; j<n; j++) { double t=LU[col*n+j]; LU[col*n+j]=LU[mr*n+j]; LU[mr*n+j]=t; }
            sign=-sign;
        }
        det*=LU[col*n+col];
        for (int r=col+1; r<n; r++) {
            double f=LU[r*n+col]/LU[col*n+col];
            for (int j=col+1; j<n; j++) LU[r*n+j] -= f*LU[col*n+j];
        }
    }
    free(LU);
    return sign*det;
}

double mat_trace(int n, const double *A)
{ if (!A||n<=0) return 0.0; double t=0.0; for(int i=0; i<n; i++) t+=A[i*n+i]; return t; }

double mat_frobenius_norm(int m, int n, const double *A)
{ if (!A||m<=0||n<=0) return 0.0; double s=0.0; for(int i=0; i<m*n; i++) s+=A[i]*A[i]; return sqrt(s); }

void mat_set_identity(int n, double *A)
{ if (!A||n<=0) return; for(int i=0; i<n*n; i++) A[i]=0.0; for(int i=0; i<n; i++) A[i*n+i]=1.0; }

void mat_scale(int m, int n, double *A, double alpha)
{ if (!A) return; for(int i=0; i<m*n; i++) A[i]*=alpha; }

void mat_copy(int m, int n, const double *src, double *dst)
{ if (!src||!dst) return; memcpy(dst,src,(size_t)m*n*sizeof(double)); }

void mat_zero(int m, int n, double *A)
{ if (!A) return; memset(A,0,(size_t)m*n*sizeof(double)); }

int mat_eigenvalues_qr(int n, double *A, double *er, double *ei, int mi, double tol)
{
    if (!A||!er||!ei||n<=0) return 0;
    (void)mi; (void)tol;
    if (n==1) { er[0]=A[0]; ei[0]=0.0; return 1; }
    double *Ac=(double*)malloc((size_t)n*n*sizeof(double));
    if (!Ac) return 0;
    memcpy(Ac,A,(size_t)n*n*sizeof(double));
    int found=0;
    for (int k=0; k<n; k++) {
        double *v=(double*)calloc((size_t)n,sizeof(double));
        double *w=(double*)calloc((size_t)n,sizeof(double));
        if (!v||!w) { free(v); free(w); break; }
        for (int i=0; i<n; i++) v[i]=(double)((i+k)%n+1);
        double lam=0.0, lam_old;
        for (int it=0; it<200; it++) {
            double norm=0.0;
            for (int i=0; i<n; i++) {
                w[i]=0.0;
                for (int j=0; j<n; j++) w[i] += Ac[i*n+j]*v[j];
                norm += w[i]*w[i];
            }
            norm=sqrt(norm); if(norm<1e-15) break;
            for (int i=0; i<n; i++) v[i]=w[i]/norm;
            double rq=0.0;
            for (int i=0; i<n; i++) {
                double av=0.0;
                for (int j=0; j<n; j++) av+=Ac[i*n+j]*v[j];
                rq+=v[i]*av;
            }
            lam_old=lam; lam=rq;
            if(fabs(lam-lam_old)<1e-10) break;
        }
        er[found]=lam; ei[found]=0.0;
        for (int i=0; i<n; i++) Ac[i*n+i]-=lam;
        found++;
        free(v); free(w);
    }
    free(Ac);
    return found;
}

int mat_power_iteration(int n, const double *A, double *x, double *lambda, int mi, double tol)
{
    if (!A||!x||!lambda||n<=0) return 0;
    double *y=(double*)malloc((size_t)n*sizeof(double));
    if (!y) return 0;
    double lam=0.0, lam_old; int iter;
    for (iter=0; iter<mi; iter++) {
        double norm=0.0;
        for (int i=0; i<n; i++) {
            y[i]=0.0;
            for (int j=0; j<n; j++) y[i]+=A[i*n+j]*x[j];
            norm+=y[i]*y[i];
        }
        norm=sqrt(norm); if(norm<1e-15) break;
        for (int i=0; i<n; i++) x[i]=y[i]/norm;
        lam_old=lam; lam=0.0;
        for (int i=0; i<n; i++) {
            double ax=0.0;
            for (int j=0; j<n; j++) ax+=A[i*n+j]*x[j];
            lam+=x[i]*ax;
        }
        if(fabs(lam-lam_old)<tol) break;
    }
    *lambda=lam; free(y);
    return iter+1;
}

int mat_lyapunov_solve(int n, const double *A, const double *Q, double *X)
{
    if (!A||!Q||!X||n<=0||n>50) return -1;
    int n2=n*n;
    double *K=(double*)calloc((size_t)n2*n2,sizeof(double));
    double *b=(double*)calloc((size_t)n2,sizeof(double));
    if (!K||!b) { free(K); free(b); return -1; }
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++) {
            int row=i*n+j;
            for (int p=0; p<n; p++) K[row*n2+(i*n+p)] += A[j*n+p];
            for (int p=0; p<n; p++) K[row*n2+(p*n+j)] += A[i*n+p];
            b[row] = -Q[i*n+j];
        }
    int ret=mat_solve(n2,K,b);
    if(ret==0) for(int i=0; i<n; i++) for(int j=0; j<n; j++) X[i*n+j]=b[i*n+j];
    free(K); free(b);
    return ret;
}

int mat_dlyap_solve(int n, const double *A, const double *Q, double *X)
{
    if (!A||!Q||!X||n<=0||n>50) return -1;
    int n2=n*n;
    double *K=(double*)calloc((size_t)n2*n2,sizeof(double));
    double *b=(double*)calloc((size_t)n2,sizeof(double));
    if (!K||!b) { free(K); free(b); return -1; }
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++) {
            int row=i*n+j;
            for (int p=0; p<n; p++)
                for (int q=0; q<n; q++)
                    K[row*n2+(p*n+q)] += A[i*n+p]*A[j*n+q];
            K[row*n2+row] -= 1.0;
            b[row] = -Q[i*n+j];
        }
    int ret=mat_solve(n2,K,b);
    if(ret==0) for(int i=0; i<n; i++) for(int j=0; j<n; j++) X[i*n+j]=b[i*n+j];
    free(K); free(b);
    return ret;
}

int mat_riccati_solve(int n, int m, const double *A, const double *B,
                       const double *Q, const double *R, double *P)
{
    if (!A||!B||!Q||!R||!P||n<=0||m<=0) return -1;
    for (int i=0; i<n*n; i++) P[i]=0.0;
    double invR=1.0/R[0];
    for (int iter=0; iter<50; iter++) {
        double *F=(double*)calloc((size_t)n*n,sizeof(double));
        if (!F) return -1;
        for (int i=0; i<n; i++)
            for (int j=0; j<n; j++) {
                double atp=0.0, pa=0.0;
                for (int k=0; k<n; k++) { atp+=A[k*n+i]*P[k*n+j]; pa+=P[i*n+k]*A[k*n+j]; }
                double pbp=0.0;
                for (int k=0; k<n; k++)
                    for (int l=0; l<n; l++) {
                        double pbr=0.0;
                        for (int r=0; r<m; r++) pbr += B[k*m+r]*B[l*m+r];
                        pbp += P[i*n+k]*pbr*invR*P[l*n+j];
                    }
                F[i*n+j] = atp + pa - 0.5*pbp + Q[i*n+j];
            }
        double nF=0.0;
        for (int i=0; i<n*n; i++) { P[i]-=0.3*F[i]; nF+=F[i]*F[i]; }
        free(F);
        if(sqrt(nF)<1e-8) return 0;
    }
    return -1;
}

int mat_transfer_func(int n, int m, int p, const double *A, const double *B,
                       const double *C, const double *D, double sigma, double omega,
                       double *G_real, double *G_imag)
{
    if (!A||!B||!C||!G_real||!G_imag||n<=0||m<=0||p<=0) return -1;
    int n2=n*n, n2x=2*n;
    double *Mreal=(double*)malloc((size_t)n2*sizeof(double));
    double *Mimag=(double*)malloc((size_t)n2*sizeof(double));
    double *Mbig=(double*)calloc((size_t)n2x*n2x,sizeof(double));
    if (!Mreal||!Mimag||!Mbig) { free(Mreal); free(Mimag); free(Mbig); return -1; }
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            Mreal[i*n+j]=-A[i*n+j]; Mimag[i*n+j]=0.0;
        }
        Mreal[i*n+i]+=sigma; Mimag[i*n+i]+=omega;
    }
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            Mbig[i*n2x+j]=Mreal[i*n+j]; Mbig[i*n2x+(n+j)]=-Mimag[i*n+j];
            Mbig[(n+i)*n2x+j]=Mimag[i*n+j]; Mbig[(n+i)*n2x+(n+j)]=Mreal[i*n+j];
        }
    }
    for (int k=0; k<m; k++) {
        double *col=(double*)calloc((size_t)n2x,sizeof(double));
        double *Mbcopy=(double*)malloc((size_t)n2x*n2x*sizeof(double));
        if (!col||!Mbcopy) { free(col); free(Mbcopy); continue; }
        memcpy(Mbcopy,Mbig,(size_t)n2x*n2x*sizeof(double));
        for (int i=0; i<n; i++) col[i]=B[i*m+k];
        if(mat_solve(n2x,Mbcopy,col)==0) {
            for (int o=0; o<p; o++) {
                double cr=0.0, ci=0.0;
                for (int j=0; j<n; j++) { cr+=C[o*n+j]*col[j]; ci+=C[o*n+j]*col[n+j]; }
                G_real[o*m+k]=cr; G_imag[o*m+k]=ci;
                if(D) G_real[o*m+k]+=D[o*m+k];
            }
        }
        free(col); free(Mbcopy);
    }
    free(Mreal); free(Mimag); free(Mbig);
    return 0;
}

int mat_freqresp_mag(int n, int m, int p, const double *A, const double *B,
                      const double *C, const double *D, double omega, double *mag)
{
    if (!mag) return -1;
    double *Gr=(double*)calloc((size_t)p*m,sizeof(double));
    double *Gi=(double*)calloc((size_t)p*m,sizeof(double));
    if (!Gr||!Gi) { free(Gr); free(Gi); return -1; }
    int ret=mat_transfer_func(n,m,p,A,B,C,D,0.0,omega,Gr,Gi);
    if(ret==0) for(int i=0; i<p*m; i++) mag[i]=sqrt(Gr[i]*Gr[i]+Gi[i]*Gi[i]);
    free(Gr); free(Gi);
    return ret;
}

double mat_sigma_max(int n, int m, int p, const double *A, const double *B,
                      const double *C, const double *D, double omega, int mi, double tol)
{
    (void)mi; (void)tol;
    double *mag=(double*)malloc((size_t)p*m*sizeof(double));
    if (!mag) return 0.0;
    int ret=mat_freqresp_mag(n,m,p,A,B,C,D,omega,mag);
    double sv=0.0;
    if(ret==0) for(int i=0; i<p*m; i++) if(mag[i]>sv) sv=mag[i];
    free(mag);
    return sv;
}
