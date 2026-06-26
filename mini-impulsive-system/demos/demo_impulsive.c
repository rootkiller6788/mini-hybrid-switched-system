/* demo_impulsive.c -- Interactive demo of impulsive systems.
 * Shows state evolution with and without impulsive control for
 * an unstable linear system.
 */
#include "impulsive_types.h"
#include "impulsive_flow.h"
#include "impulsive_jump.h"
#include "impulsive_lyapunov.h"
#include "impulsive_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(void) {
    printf("=== Impulsive System Demo ===\n");
    printf("Comparing uncontrolled vs impulsively controlled unstable system.\n\n");
    int n=2; double A[4]={1.0,0.0,0.0,1.0};
    ImpFlowLinear *flow=affine_flow_create(A,n);
    ImpSolverConfig cfg=imp_solver_config_default();
    cfg.h_init=0.01;
    double x[2]={1.0,1.0},t=0.0,T=3.0;
    printf("Uncontrolled (grows exponentially):\n");
    while(t<T) {
        double xn[2]; imp_solver_rk4_step((ImpVectorField)affine_flow_eval,flow,t,x,n,cfg.h_init,xn);
        t+=cfg.h_init; x[0]=xn[0]; x[1]=xn[1];
        if(fmod(t,0.5)<cfg.h_init) printf("  t=%.1f: [%.3f, %.3f] (norm=%.3f)\n",t,x[0],x[1],sqrt(x[0]*x[0]+x[1]*x[1]));
    }
    double J[4]={-0.9,0.0,0.0,-0.9};
    ImpJumpLinear *jump=imp_jump_linear_create(J,n);
    x[0]=1.0;x[1]=1.0;t=0.0;
    printf("\nImpulsively controlled (stabilized):\n");
    for(int k=0;k<10;k++) {
        double xn[2]; imp_solver_rk4_step((ImpVectorField)affine_flow_eval,flow,t,x,n,0.3,xn);
        t+=0.3; x[0]=xn[0];x[1]=xn[1];
        double xb[2]={x[0],x[1]}; imp_jump_linear_eval(t,xb,n,x,jump);
        printf("  t=%.1f: [%.3f, %.3f] (norm=%.3f)\n",t,x[0],x[1],sqrt(x[0]*x[0]+x[1]*x[1]));
    }
    affine_flow_free(flow); imp_jump_linear_free(jump);
    return 0;
}
