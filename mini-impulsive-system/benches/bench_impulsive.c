#include "impulsive_types.h"
#include "impulsive_flow.h"
#include "impulsive_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
int main(void) {
    int n=3; double A[9]={-1,0,0,0,-2,0,0,0,-3};
    ImpFlowLinear *fl=affine_flow_create(A,n);
    ImpSolverConfig cfg=imp_solver_config_default();
    cfg.h_init=0.001; cfg.method=IMP_SOLVER_RK4;
    double x[3]={1,2,3},xn[3];
    int steps=100000;
    clock_t start=clock();
    for(int i=0;i<steps;i++)
        imp_solver_rk4_step((ImpVectorField)affine_flow_eval,fl,0.0,x,n,cfg.h_init,xn);
    clock_t end=clock();
    double dt=(double)(end-start)/CLOCKS_PER_SEC;
    printf("RK4 bench: %d steps in %.4fs (%.1f steps/s)\n",steps,dt,(double)steps/dt);
    affine_flow_free(fl); return 0;
}
