/*
 * test_impulsive.c -- Unit tests for mini-impulsive-system
 * Tests all core APIs: types, flow, jump, Lyapunov, solver, sync, control.
 */
#include "impulsive_types.h"
#include "impulsive_flow.h"
#include "impulsive_jump.h"
#include "impulsive_lyapunov.h"
#include "impulsive_solver.h"
#include "impulsive_sync.h"
#include "impulsive_control.h"
#include "impulsive_analysis.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#define N 3
#define TOL 1e-6

static int tf(double t, const double *x, int n, double *dxdt, void *ctx)
{ (void)t; (void)ctx; for(int i=0;i<n;i++) dxdt[i]=x[i]; return 0; }

static int tj(double t, const double *x, int n, double *xa, void *ctx)
{ (void)t; (void)ctx; memcpy(xa,x,(size_t)n*sizeof(double)); return 0; }

static int tr=0, tp=0;
#define T(nm) do{tr++;printf("  %s ... ",nm);}while(0)
#define P() do{tp++;printf("PASS\n");}while(0)
#define F(m) do{printf("FAIL: %s\n",m);}while(0)

int main(void) {
    printf("=== mini-impulsive-system test suite ===\n\n");

    T("create_free_seq");
    { ImpTimeSeq *s=imp_time_seq_create(10,0.0,10.0); assert(s); imp_time_seq_free(s); P(); }

    T("periodic_seq");
    { ImpTimeSeq *s=imp_time_seq_create_periodic(0.0,10.0,0.5);
      assert(s&&s->count==21); assert(fabs(s->times[0])<TOL);
      assert(fabs(s->times[1]-0.5)<TOL); imp_time_seq_free(s); P(); }

    T("add_impulse_time");
    { ImpTimeSeq *s=imp_time_seq_create(4,0.0,10.0);
      assert(imp_time_seq_add(s,1.0)==0); assert(imp_time_seq_add(s,3.0)==0);
      assert(imp_time_seq_add(s,7.0)==0); assert(s->count==3);
      assert(imp_time_seq_add(s,0.5)==-1); imp_time_seq_free(s); P(); }

    T("dwell_stats");
    { ImpTimeSeq *s=imp_time_seq_create_periodic(0.0,10.0,2.0);
      ImpDwellStats st; assert(imp_time_seq_get_dwell_stats(s,&st)==0);
      assert(fabs(st.min_dwell-2.0)<TOL); assert(fabs(st.max_dwell-2.0)<TOL);
      imp_time_seq_free(s); P(); }

    T("find_index");
    { double t[]={0.0,2.0,4.0,6.0,8.0};
      ImpTimeSeq *s=imp_time_seq_create_aperiodic(t,5,0.0,10.0);
      assert(s); assert(imp_time_seq_find_index(s,1.5)==0);
      assert(imp_time_seq_find_index(s,3.0)==1);
      assert(imp_time_seq_find_index(s,9.0)==4); imp_time_seq_free(s); P(); }

    T("system_create_validate");
    { ImpTimeSeq *s=imp_time_seq_create_periodic(0.0,10.0,1.0);
      ImpSystem *sys=imp_system_create(N,tf,tj,s,NULL);
      assert(sys&&imp_system_validate(sys)); imp_system_free(sys);
      imp_time_seq_free(s); P(); }

    T("lorenz_flow");
    { ImpFlowLorenz *l=imp_flow_lorenz_create(10.0,28.0,8.0/3.0);
      assert(l); double x[3]={1,1,1},d[3];
      assert(imp_flow_lorenz_eval(0.0,x,3,d,l)==0);
      assert(fabs(d[0]-0.0)<TOL); imp_flow_lorenz_free(l); P(); }

    T("vdp_flow");
    { ImpFlowVanDerPol *v=imp_flow_vdp_create(2.0);
      assert(v); double x[2]={2,0},d[2];
      assert(imp_flow_vdp_eval(0.0,x,2,d,v)==0);
      assert(fabs(d[0])<TOL); imp_flow_vdp_free(v); P(); }

    T("linear_flow");
    { double A[4]={0,1,-1,0}; ImpFlowLinear *fl=affine_flow_create(A,2);
      assert(fl); double x[2]={1,0},d[2];
      assert(affine_flow_eval(0.0,x,2,d,fl)==0);
      assert(fabs(d[0])<TOL); assert(fabs(d[1]+1.0)<TOL);
      affine_flow_free(fl); P(); }

    T("linear_jump");
    { double J[4]={-0.5,0,0,-0.5}; ImpJumpLinear *jl=imp_jump_linear_create(J,2);
      assert(jl); double xb[2]={2,4},xa[2];
      assert(imp_jump_linear_eval(0.0,xb,2,xa,jl)==0);
      assert(fabs(xa[0]-1.0)<TOL); assert(fabs(xa[1]-2.0)<TOL);
      imp_jump_linear_free(jl); P(); }

    T("projection_jump");
    { double lb[2]={-1,-1},ub[2]={1,1};
      ImpJumpProject *jp=imp_jump_project_create(lb,ub,2);
      assert(jp); double xb[2]={2,-3},xa[2];
      assert(imp_jump_project_eval(0.0,xb,2,xa,jp)==0);
      assert(fabs(xa[0]-1.0)<TOL); assert(fabs(xa[1]+1.0)<TOL);
      imp_jump_project_free(jp); P(); }

    T("quadratic_lyapunov");
    { double P[4]={1,0,0,2}; ImpLyapQuad *ly=imp_lyap_quad_create(P,2);
      assert(ly); double x[2]={3,2};
      assert(fabs(imp_lyap_quad_V(x,2,ly)-17.0)<TOL); imp_lyap_quad_free(ly); P(); }

    T("dwell_time_bound");
    { assert(fabs(imp_lyap_dwell_time_bound(1.0,0.5)-log(2.0))<0.01); P(); }

    T("lyapunov_eq_solver");
    { double A[4]={-1,0,0,-2},Q[4]={1,0,0,1},P[4];
      assert(imp_lyap_solve_lyapunov_eq(A,Q,P,2)==0);
      assert(P[0]>0&&P[3]>0); P(); }

    T("rk4_step");
    { double A[4]={0,1,-1,0}; ImpFlowLinear *fl=affine_flow_create(A,2);
      double x[2]={1,0},xn[2];
      assert(imp_solver_rk4_step((ImpVectorField)affine_flow_eval,fl,0.0,x,2,0.01,xn)==0);
      assert(fabs(xn[0]-cos(0.01))<0.001); assert(fabs(xn[1]+sin(0.01))<0.001);
      affine_flow_free(fl); P(); }

    T("scalar_sync_jump");
    { ImpSyncJumpScalar *sjs=imp_sync_jump_scalar_create(0.8,2);
      double eb[2]={1,-2},ea[2];
      assert(imp_sync_jump_scalar_eval(0.0,eb,2,ea,sjs)==0);
      assert(fabs(ea[0]-0.2)<TOL); assert(fabs(ea[1]+0.4)<TOL);
      imp_sync_jump_scalar_free(sjs); P(); }

    T("sync_condition");
    { assert(imp_sync_check_synchronization_condition(0.5,1.0,0.8)); P(); }

    T("control_effort");
    { double u[3]={1,2,2}; assert(fabs(imp_control_effort_L2(u,3)-3.0)<TOL); P(); }

    T("dwell_bound_check");
    { double t[]={0,5,12,20}; ImpTimeSeq *s=imp_time_seq_create_aperiodic(t,4,0,20);
      assert(imp_analysis_dwell_satisfies_bound(s,3.0));
      assert(!imp_analysis_dwell_satisfies_bound(s,10.0));
      imp_time_seq_free(s); P(); }

    T("energy_analysis");
    { double x[3]={3,4,0}; assert(fabs(imp_analysis_energy(x,3)-12.5)<TOL); P(); }

    T("invariant_set");
    { ImpInvariantSet *inv=imp_analysis_invariant_set_create(2);
      double x[2]={0.5,0.5}; assert(imp_analysis_is_in_set(inv,x));
      imp_analysis_invariant_set_free(inv); P(); }

    T("monodromy");
    { double A[4]={-1,0,0,-2}; ImpFlowLinear *fl=affine_flow_create(A,2);
      ImpTimeSeq *s=imp_time_seq_create_periodic(0.0,1.0,1.0);
      ImpSystem *sys=imp_system_create(2,(ImpVectorField)affine_flow_eval,tj,s,fl);
      double Phi[4];
      ImpSolverConfig cfg=imp_solver_config_default();
      assert(imp_analysis_monodromy_matrix(sys,&cfg,Phi,2)==0);
      imp_system_free(sys); imp_time_seq_free(s); affine_flow_free(fl);
      P(); }

    printf("\n=== Results: %d/%d passed ===\n",tp,tr);
    return (tp==tr)?0:1;
}
