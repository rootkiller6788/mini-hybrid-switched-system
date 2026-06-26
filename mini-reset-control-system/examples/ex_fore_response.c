/* ex_fore_response.c - First-Order Reset Element Response Demo
 *
 * Demonstrates the FORE response characteristics at different reset
 * ratios rho. As rho varies from 0 (full reset, Clegg-like) to 0.9
 * (nearly linear), the phase advantage of reset control diminishes
 * linearly.
 *
 * This example shows:
 *   - FORE step response for rho = 0.0, 0.3, 0.6, 0.9
 *   - Effect of reset ratio on steady-state behavior
 *   - Reset event counting per configuration
 *
 * Ref: Horowitz & Rosenbaum (1975); Banos & Barreiro (2012) Section 3.3
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "reset_element.h"

void simulate_fore(double K, double tau, double rho, double dt, double t_end)
{
    ForeElement *fe = fore_create(K, tau, rho);
    if (!fe) {
        printf("  ERROR: Failed to create FORE with rho=%.2f\n", rho);

        return;
    }

    int n_steps = (int)(t_end / dt);
    double y = 0.0, e_prev = 0.0;
    double y_final = 0.0;
    int reset_count = 0;

    /* Apply step input of 1.0 */
    for (int k = 0; k < n_steps; k++) {
        double e = 1.0;
        y = fore_step(fe, dt, e, e_prev);
        e_prev = e;
        if (k == n_steps - 1) y_final = y;
    }

    /* Count resets via the base system */
    const ResetInterval *ri = reset_get_interval_stats(fe->base);
    if (ri) reset_count = ri->n_resets;

    printf("  rho=%.2f | steady-state output=%.4f | resets=%d\n",

rho, y_final, reset_count);


    fore_free(fe);
}

int main(void)
{
    printf("=== FORE (First-Order Reset Element) Response Demo ===\n");


    double K = 1.0;      /* DC gain */
    double tau = 0.5;     /* time constant */
    double dt = 0.001;
    double t_end = 5.0;

    printf("Parameters: K=%.1f, tau=%.2f, dt=%.4f, T_end=%.1f\n",

K, tau, dt, t_end);

    printf("Step response for different reset ratios:\n");


    simulate_fore(K, tau, 0.0, dt, t_end);  /* full reset (Clegg-like) */
    simulate_fore(K, tau, 0.3, dt, t_end);
    simulate_fore(K, tau, 0.6, dt, t_end);
    simulate_fore(K, tau, 0.9, dt, t_end);  /* nearly linear */

    printf("Analysis:\n");

    printf("  - rho=0: Full reset, maximum phase advantage\n");

    printf("  - rho=0.9: Nearly linear, minimal reset effect\n");

    printf("  - Higher rho = less phase lead but smoother response\n");


    printf("=== Demo Complete ===\n");

    return 0;
}
