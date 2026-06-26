/* test_elements.c - Tests for reset control elements */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "reset_element.h"

static int passed=0, failed=0;
#define T(name) printf("  TEST: %s ... ",name)
#define P() do{printf("PASS\n");passed++;}while(0)
#define F(m) do{printf("FAIL: %s\n",m);failed++;}while(0)
#define C(c) do{if(!(c)){F(#c);return;}}while(0)

void test_clegg_step(void) {
    T("clegg_step");
    CleggIntegrator *ci=clegg_create();
    C(ci!=NULL);
    double y=clegg_step(ci,0.01,1.0,0.5);
    C(fabs(y-0.01)<1e-5);
    y=clegg_step(ci,0.01,-0.1,1.0);
    C(fabs(y)<1e-2);
    clegg_free(ci); P();
}

void test_fore_step(void) {
    T("fore_step");
    ForeElement *fe=fore_create(1.0,0.1,0.5);
    C(fe!=NULL);
    double y=fore_step(fe,0.01,1.0,0.5);
    C(fabs(y-0.0)<0.2);
    y=fore_step(fe,0.01,-0.5,1.0);
    C(fabs(y)<0.2);
    fore_free(fe); P();
}

void test_reset_pid(void) {
    T("reset_pid");
    ResetPID *pid=reset_pid_create(1.0,0.5,0.1,0.01,0.0);
    C(pid!=NULL);
    double y=reset_pid_step(pid,0.01,1.0,0.5);
    C(fabs(y)>0.5);
    double yp,yi,yd;
    reset_pid_get_components(pid,&yp,&yi,&yd);
    C(fabs(yp-1.0)<0.1);
    reset_pid_free(pid); P();
}

void test_sore_step(void) {
    T("sore_step");
    SoreElement *se=sore_create(10.0,0.7,0.0,0.0);
    C(se!=NULL);
    double y=sore_step(se,0.001,1.0,0.5);
    C(isfinite(y));
    sore_free(se); P();
}

void test_clegg_df(void) {
    T("clegg_df");
    double mag=clegg_df_magnitude(1.0);
    C(fabs(mag-1.62)<0.01);
    double ph=clegg_df_phase();
    C(fabs(ph-(-0.905))<0.01);
    P();
}

void test_fore_df(void) {
    T("fore_df");
    ForeElement *fe=fore_create(1.0,0.5,0.2);
    C(fe!=NULL);
    double mag=fore_df_magnitude(fe,1.0);
    C(mag>0.0);
    double ph=fore_df_phase(fe,1.0);
    C(ph>-2.0 && ph<2.0);
    fore_free(fe); P();
}

void test_leadlag(void) {
    T("leadlag");
    ResetLeadLag *rll=reset_leadlag_create(1.0,0.5,0.1,0.3,true);
    C(rll!=NULL);
    double y=reset_leadlag_step(rll,0.01,1.0,0.5);
    C(isfinite(y));
    double mag,ph;
    reset_leadlag_freqresp(rll,1.0,&mag,&ph);
    C(mag>0.0);
    reset_leadlag_free(rll); P();
}

int main(void) {
    printf("=== Reset Element Tests ===\n\n");
    test_clegg_step();
    test_fore_step();
    test_reset_pid();
    test_sore_step();
    test_clegg_df();
    test_fore_df();
    test_leadlag();
    printf("\n=== Results: %d passed, %d failed ===\n",passed,failed);
    return failed>0?1:0;
}
