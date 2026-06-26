/**
 * @file example_cruise_control.c
 * @brief L6 KP4: Adaptive Cruise Control — Automotive Hybrid System
 *
 * Demonstrates a three-mode hybrid system for vehicle cruise control:
 *   CRUISE:  maintain set speed
 *   FOLLOW:  maintain safe time headway behind lead vehicle
 *   BRAKE:   emergency deceleration to avoid collision
 *
 * Mode switching is state-dependent based on the gap to the lead vehicle.
 * This is a practical example of hybrid systems in automotive control.
 */
#include "hss_simulation.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Adaptive Cruise Control — Automotive Hybrid System (L6) ===\n\n");

    /* Lead vehicle ahead at 50m, moving at 28 m/s (~100 km/h) */
    /* Ego vehicle at 30 m/s (~108 km/h) with set speed 35 m/s (~126 km/h) */
    HSS_CruiseControl cc = hss_cruise_control_init(
        30.0,    /* ego_vel: 30 m/s */
        50.0,    /* lead_pos: 50 m ahead */
        28.0,    /* lead_vel: 28 m/s */
        35.0,    /* v_set: 35 m/s desired */
        2.0,     /* headway: 2 seconds */
        3.0      /* min_gap: 3 m at standstill */
    );

    double desired_gap = cc.min_gap + cc.time_headway * cc.ego_velocity;
    printf("Parameters:\n");
    printf("  Set speed: %.1f m/s (%.0f km/h)\n", cc.v_set, cc.v_set * 3.6);
    printf("  Headway: %.1f s, Min gap: %.1f m\n",
           cc.time_headway, cc.min_gap);
    printf("  Initial gap: %.1f m (desired: %.1f m)\n",
           cc.lead_position - cc.ego_position, desired_gap);
    printf("  Initial mode: %s\n",
           cc.mode == 0 ? "CRUISE" : cc.mode == 1 ? "FOLLOW" : "BRAKE");
    printf("\n");

    double dt = 0.05;
    double sim_time = 30.0;
    double print_dt = 2.0;
    double next_print = 0.0;

    printf("Simulation trace:\n");
    printf("  %8s  %8s  %8s  %8s  %8s  %8s\n",
           "Time", "EgoV", "LeadV", "Gap", "Mode", "CollRisk");

    while (cc.time < sim_time) {
        int mode = hss_cruise_control_step(&cc, dt);

        if (cc.time >= next_print) {
            double gap = cc.lead_position - cc.ego_position;
            const char *mode_str = (mode == 0) ? "CRUISE" :
                                   (mode == 1) ? "FOLLOW" : "BRAKE ";
            printf("  %8.2f  %8.2f  %8.2f  %8.2f  %8s  %8s\n",
                   cc.time, cc.ego_velocity, cc.lead_velocity,
                   gap, mode_str, cc.collision_risk ? "YES" : "no");
            next_print += print_dt;
        }
    }

    printf("\nFinal state:\n");
    printf("  Ego position: %.1f m, velocity: %.2f m/s\n",
           cc.ego_position, cc.ego_velocity);
    printf("  Lead position: %.1f m\n", cc.lead_position);
    printf("  Gap: %.1f m (desired: %.1f m)\n",
           cc.lead_position - cc.ego_position,
           cc.min_gap + cc.time_headway * cc.ego_velocity);
    printf("  Gap error: %.2f m\n", cc.gap_error);
    printf("  Final mode: %s\n",
           cc.mode == 0 ? "CRUISE" :
           cc.mode == 1 ? "FOLLOW" : "BRAKE");
    printf("  Collision risk: %s\n",
           cc.collision_risk ? "YES" : "No");

    return 0;
}
