/* reset_simulation.c - Hybrid simulation engine with ZC detection */
#include "reset_simulation.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <float.h>

SimConfig sim_config_default(void)
{
    SimConfig cfg;
    cfg.method = SIM_METHOD_RK4;
    cfg.dt_min = 1e-8;
    cfg.dt_max = 0.01;
    cfg.dt_fixed = 0.001;
    cfg.t_start = 0.0;
    cfg.t_end = 10.0;
    cfg.zc_tol = 1e-8;
    cfg.max_resets = 10000;
    cfg.max_steps = 1000000;
    cfg.zc_refine_tol = 1e-10;
    cfg.zc_max_iter = 20;
    cfg.adaptive = false;
    cfg.detect_zeno = true;
    cfg.zeno_threshold = 1000.0;
    return cfg;
}

SimState* sim_state_create(ResetSystem *rsys)
{
    if (!rsys) return NULL;
    SimState *st=(SimState*)calloc(1,sizeof(SimState));
    if (!st) return NULL;
    st->rsys=rsys;
    st->t=0.0;
    st->x=(double*)calloc((size_t)rsys->nc,sizeof(double));
    if (!st->x) { free(st); return NULL; }
    memcpy(st->x,rsys->x_c,(size_t)rsys->nc*sizeof(double));
    st->step_count=0;
    st->reset_count=0;
    st->t_last_reset=-1.0;
    st->dt_current=0.001;
    st->in_flow=true;
    st->zeno_warning=false;
    st->log_capacity=1000;
    st->log_count=0;
    st->reset_log=(ResetMemory*)calloc((size_t)st->log_capacity,sizeof(ResetMemory));
    if (!st->reset_log) { free(st->x); free(st); return NULL; }
    return st;
}

void sim_state_free(SimState *st)
{
    if (!st) return;
    free(st->x);
    if (st->reset_log) {
        for (int i=0; i<st->log_count; i++) {
            free(st->reset_log[i].x_before);
            free(st->reset_log[i].x_after);
        }
        free(st->reset_log);
    }
    if (st->u_vec) free(st->u_vec);
    free(st);
}

void sim_set_initial_state(SimState *st, const double *x0)
{
    if (!st||!x0) return;
    memcpy(st->x,x0,(size_t)st->rsys->nc*sizeof(double));
    memcpy(st->rsys->x_c,x0,(size_t)st->rsys->nc*sizeof(double));
}

double sim_integrate_flow(SimState *st, const SimConfig *cfg, double u)
{
    if (!st||!cfg||!st->rsys||!st->rsys->flow) return 0.0;
    st->u=u;
    double dt=cfg->adaptive?st->dt_current:cfg->dt_fixed;
    int n=st->rsys->nc;
    ResetLinearBase *fl=st->rsys->flow;
    double *x0=(double*)malloc((size_t)n*sizeof(double));
    double *k1=(double*)malloc((size_t)n*sizeof(double));
    double *k2=(double*)malloc((size_t)n*sizeof(double));
    double *k3=(double*)malloc((size_t)n*sizeof(double));
    double *k4=(double*)malloc((size_t)n*sizeof(double));
    double *xt=(double*)malloc((size_t)n*sizeof(double));
    if (!x0||!k1||!k2||!k3||!k4||!xt) {
        free(x0); free(k1); free(k2); free(k3); free(k4); free(xt);
        return dt;
    }
    memcpy(x0,st->x,(size_t)n*sizeof(double));
    double uv[1]; uv[0]=u;
    if (cfg->method==SIM_METHOD_EULER) {
        reset_eval_flow_deriv(fl,x0,uv,k1);
        for (int i=0; i<n; i++) st->x[i]=x0[i]+dt*k1[i];
    } else if (cfg->method==SIM_METHOD_HEUN) {
        reset_eval_flow_deriv(fl,x0,uv,k1);
        for (int i=0; i<n; i++) xt[i]=x0[i]+dt*k1[i];
        reset_eval_flow_deriv(fl,xt,uv,k2);
        for (int i=0; i<n; i++) st->x[i]=x0[i]+0.5*dt*(k1[i]+k2[i]);
    } else {
        /* RK4 */
        reset_eval_flow_deriv(fl,x0,uv,k1);
        for (int i=0; i<n; i++) xt[i]=x0[i]+0.5*dt*k1[i];
        reset_eval_flow_deriv(fl,xt,uv,k2);
        for (int i=0; i<n; i++) xt[i]=x0[i]+0.5*dt*k2[i];
        reset_eval_flow_deriv(fl,xt,uv,k3);
        for (int i=0; i<n; i++) xt[i]=x0[i]+dt*k3[i];
        reset_eval_flow_deriv(fl,xt,uv,k4);
        for (int i=0; i<n; i++) st->x[i]=x0[i]+(dt/6.0)*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);
    }
    st->t+=dt;
    st->step_count++;
    memcpy(st->rsys->x_c,st->x,(size_t)n*sizeof(double));
    st->rsys->t=st->t;
    free(x0); free(k1); free(k2); free(k3); free(k4); free(xt);
    return dt;
}

bool sim_check_reset(SimState *st, const SimConfig *cfg, double e_now, double e_prev)
{
    if (!st||!cfg||!st->rsys) return false;
    ResetResult rr=reset_check_and_apply(st->rsys,e_now,e_prev);
    if (rr==RESET_OK) {
        st->reset_count++;
        st->t_last_reset=st->t;
        /* Log the reset */
        if (st->log_count<st->log_capacity) {
            int idx=st->log_count;
            st->reset_log[idx].x_before=(double*)malloc((size_t)st->rsys->nc*sizeof(double));
            st->reset_log[idx].x_after=(double*)malloc((size_t)st->rsys->nc*sizeof(double));
            if (st->reset_log[idx].x_before&&st->reset_log[idx].x_after) {
                memcpy(st->reset_log[idx].x_before,st->x,(size_t)st->rsys->nc*sizeof(double));
                memcpy(st->reset_log[idx].x_after,st->rsys->x_c,(size_t)st->rsys->nc*sizeof(double));
                st->reset_log[idx].t_reset=st->t;
                st->reset_log[idx].e_at_reset=e_now;
                st->reset_log[idx].reset_seq=st->reset_count;
                st->log_count++;
            }
        }
        memcpy(st->x,st->rsys->x_c,(size_t)st->rsys->nc*sizeof(double));
        /* Check Zeno */
        if (cfg->detect_zeno) {
            double interval=st->t-st->t_last_reset;
            if (st->reset_count>10&&interval<cfg->zc_tol) st->zeno_warning=true;
        }
        return true;
    }
    return false;
}

SimResult* sim_run(SimState *st, const SimConfig *cfg, const double *u_vals, int n_u)
{
    if (!st||!cfg||!u_vals||n_u<=0) return NULL;
    SimResult *res=(SimResult*)calloc(1,sizeof(SimResult));
    if (!res) return NULL;
    int cap=n_u+10;
    res->t_history=(double*)malloc((size_t)cap*sizeof(double));
    res->x_history=(double*)malloc((size_t)cap*st->rsys->nc*sizeof(double));
    res->y_history=(double*)malloc((size_t)cap*st->rsys->flow->p*sizeof(double));
    if (!res->t_history||!res->x_history||!res->y_history) { sim_result_free(res); return NULL; }
    int hist=0; res->t_history[0]=st->t;
    for (int i=0; i<st->rsys->nc; i++) res->x_history[0*st->rsys->nc+i]=st->x[i];
    if (st->rsys->flow->p==1) res->y_history[0]=st->rsys->flow->C[0]*st->x[0];
    hist=1;
    for (int k=1; k<n_u&&k<cap; k++) {
        double dt_hint=0.0;
        int sub_steps=1;
        double dt=(cfg->t_end-cfg->t_start)/(double)n_u;
        sub_steps=(dt>cfg->dt_fixed)?(int)(dt/cfg->dt_fixed)+1:1;
        double sub_dt=dt/(double)sub_steps;
        for (int sub=0; sub<sub_steps; sub++) {
            double dt_used=sim_integrate_flow(st,cfg,u_vals[k]);
            (void)dt_used;
            dt_hint+=sub_dt;
        }
        /* Check reset */
        double e_prev=st->rsys->e_prev;
        double e_now=u_vals[k];
        sim_check_reset(st,cfg,e_now,e_prev);
        /* Record */
        res->t_history[hist]=st->t;
        for (int i=0; i<st->rsys->nc; i++) res->x_history[hist*st->rsys->nc+i]=st->x[i];
        if (st->rsys->flow->p==1) res->y_history[hist]=st->rsys->flow->C[0]*st->x[0];
        hist++;
    }
    res->n_steps=st->step_count;
    res->n_resets=st->reset_count;
    res->t_final=st->t;
    res->completed=true;
    res->zeno_detected=st->zeno_warning;
    res->history_len=hist;
    return res;
}

SimResult* sim_run_feedback(SimState *st, const SimConfig *cfg, const double *r_vals, int n_r)
{
    if (!st||!cfg) return NULL;
    SimResult *res=(SimResult*)calloc(1,sizeof(SimResult));
    if (!res) return NULL;
    int cap=n_r+10;
    res->t_history=(double*)malloc((size_t)cap*sizeof(double));
    res->x_history=(double*)malloc((size_t)cap*st->rsys->nc*sizeof(double));
    res->y_history=(double*)malloc((size_t)cap*sizeof(double));
    if (!res->t_history||!res->x_history||!res->y_history) { sim_result_free(res); return NULL; }
    int hist=0;
    res->t_history[0]=st->t;
    res->y_history[0]=0.0;
    hist=1;
    for (int k=1; k<n_r&&k<cap; k++) {
        double dt=(cfg->t_end-cfg->t_start)/(double)n_r;
        int sub=(int)(dt/cfg->dt_fixed)+1;
        double sd=dt/(double)sub; (void)sd;
        for (int ss=0; ss<sub; ss++) {
            double u_err=r_vals[k]-res->y_history[hist-1];
            double dt_used=sim_integrate_flow(st,cfg,u_err);
            (void)dt_used;
            double e_prev=st->rsys->e_prev;
            sim_check_reset(st,cfg,u_err,e_prev);
            st->rsys->e_prev=u_err;
        }
        res->t_history[hist]=st->t;
        for (int i=0; i<st->rsys->nc; i++) res->x_history[hist*st->rsys->nc+i]=st->x[i];
        res->y_history[hist]=st->rsys->flow->C[0]*st->x[0];
        hist++;
    }
    res->n_steps=st->step_count;
    res->n_resets=st->reset_count;
    res->t_final=st->t;
    res->completed=true;
    res->zeno_detected=st->zeno_warning;
    res->history_len=hist;
    return res;
}

void sim_result_free(SimResult *r)
{
    if (!r) return;
    free(r->t_history); free(r->x_history); free(r->y_history);
    free(r);
}

bool sim_detect_zc(double e0, double e1, double t0, double t1, double *zc_time)
{
    if (!zc_time) return false;
    if (e0*e1>0.0) return false;
    if (fabs(e0)<1e-15) { *zc_time=t0; return true; }
    if (fabs(e1)<1e-15) { *zc_time=t1; return true; }
    *zc_time=t0-e0*(t1-t0)/(e1-e0);
    return true;
}

double sim_refine_zc(sim_error_func f_eval, void *ctx, double ta, double tb,
                      double tol, int max_iter)
{
    double fa=f_eval(ta,ctx), fb=f_eval(tb,ctx);
    if (fa*fb>0.0) return ta;
    double tc=ta;
    for (int i=0; i<max_iter; i++) {
        tc=(ta+tb)/2.0;
        double fc=f_eval(tc,ctx);
        if (fabs(fc)<tol||(tb-ta)/2.0<tol) break;
        if (fa*fc<0.0) { tb=tc; fb=fc; } else { ta=tc; fa=fc; }
    }
    return tc;
}

double sim_detect_zeno(const SimState *st, double wsize, double zthresh)
{
    if (!st||st->reset_count<10) return INFINITY;
    (void)wsize; (void)zthresh;
    if (st->zeno_warning) return st->t;
    return INFINITY;
}

void sim_print_history(const SimResult *r, FILE *fp)
{
    if (!r||!fp) return;
    for (int i=0; i<r->history_len; i++)
        fprintf(fp,"%g,%g,%g\n",r->t_history[i],r->y_history[i],r->x_history[i]);
}

const ResetMemory* sim_get_reset_log(const SimState *st, int *count)
{
    if (!st) return NULL;
    if (count) *count=st->log_count;
    return st->reset_log;
}

double sim_compute_rms_error(const SimResult *r, const double *ref, int n_ref)
{
    if (!r||!ref||n_ref<=0) return 0.0;
    double sum=0.0; int cnt=0;
    for (int i=0; i<r->history_len&&i<n_ref; i++) {
        double err=r->y_history[i]-ref[i];
        sum+=err*err; cnt++;
    }
    return cnt>0?sqrt(sum/(double)cnt):0.0;
}

double sim_settling_time(const SimResult *r, double tol_pct)
{
    if (!r||r->history_len<2) return INFINITY;
    double y_final=r->y_history[r->history_len-1];
    double tol=fabs(y_final)*tol_pct/100.0;
    if (tol<1e-10) tol=0.01;
    for (int i=r->history_len-1; i>=0; i--) {
        if (fabs(r->y_history[i]-y_final)>tol) {
            return (i+1<r->history_len)?r->t_history[i+1]:r->t_history[r->history_len-1];
        }
    }
    return 0.0;
}
