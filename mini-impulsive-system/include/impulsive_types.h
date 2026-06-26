#ifndef IMPULSIVE_TYPES_H
#define IMPULSIVE_TYPES_H
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#define IMP_PI         3.14159265358979323846
#define IMP_EPS        1e-12
#define IMP_SQRT_EPS   1e-6
#define IMP_MAX_DIM    32
#define IMP_MAX_IMPULSES 4096
#define IMP_MAX_ITER   10000
#define IMP_HUGE       1e308
typedef enum { IMP_TRIGGER_TIME=0, IMP_TRIGGER_STATE=1, IMP_TRIGGER_MIXED=2, IMP_TRIGGER_EVENT=3 } ImpTriggerType;
typedef enum { IMP_JUMP_LINEAR=0, IMP_JUMP_AFFINE=1, IMP_JUMP_NONLINEAR=2, IMP_JUMP_PROJECTION=3, IMP_JUMP_HARD_RESET=4, IMP_JUMP_IMPULSE_MAP=5 } ImpJumpType;
typedef enum { IMP_STABLE=0, IMP_ASYMPTOTIC=1, IMP_EXPONENTIAL=2, IMP_PRACTICAL=3, IMP_FINITE_TIME=4, IMP_INPUT_TO_STATE=5 } ImpStabilityType;
typedef struct ImpSystem ImpSystem;
typedef struct ImpSolution ImpSolution;
typedef struct ImpLyapunovFn ImpLyapunovFn;
typedef int (*ImpVectorField)(double t, const double *x, int n, double *dxdt, void *ctx);
typedef int (*ImpJumpMap)(double t_k, const double *x_before, int n, double *x_after, void *ctx);
typedef bool (*ImpGuardFunction)(double t, const double *x, int n, void *ctx);
typedef double (*ImpEventFunction)(double t, const double *x, int n, void *ctx);
typedef struct { double *times; int count; int capacity; double t0; double T; bool is_periodic; double period; } ImpTimeSeq;
struct ImpSystem { int n; ImpVectorField f; ImpJumpMap I; ImpTimeSeq *seq; ImpTriggerType trigger; ImpGuardFunction guard; ImpEventFunction event; void *ctx; char name[64]; };
struct ImpSolution { double *t; double *x; double *x_jump_before; double *x_jump_after; int npts; int npts_cap; int n; int nimp; int nimp_cap; double t0; double T; };
typedef struct { double min_dwell; double max_dwell; double avg_dwell; double total_time; int num_dwells; } ImpDwellStats;
typedef struct { double (*g)(double t, double u, void *ctx); double (*psi)(double u, void *ctx); void *ctx; } ImpComparisonSys;

ImpTimeSeq* imp_time_seq_create(int capacity, double t0, double T);
void imp_time_seq_free(ImpTimeSeq *seq);
int imp_time_seq_add(ImpTimeSeq *seq, double tau);
int imp_time_seq_complete(ImpTimeSeq *seq);
ImpTimeSeq* imp_time_seq_create_periodic(double t0, double T, double period);
ImpTimeSeq* imp_time_seq_create_aperiodic(const double *times, int count, double t0, double T);
int imp_time_seq_get_dwell_stats(const ImpTimeSeq *seq, ImpDwellStats *stats);
int imp_time_seq_find_index(const ImpTimeSeq *seq, double t);
ImpTimeSeq* imp_time_seq_clone(const ImpTimeSeq *src);
ImpSystem* imp_system_create(int n, ImpVectorField f, ImpJumpMap I, ImpTimeSeq *seq, void *ctx);
void imp_system_free(ImpSystem *sys);
bool imp_system_validate(const ImpSystem *sys);
ImpSolution* imp_solution_create(int npts_cap, int nimp_cap, int n);
void imp_solution_free(ImpSolution *sol);
int imp_solution_resize(ImpSolution *sol);
int imp_solution_add_point(ImpSolution *sol, double t, const double *x);
int imp_solution_add_jump(ImpSolution *sol, const double *x_before, const double *x_after);
#endif
