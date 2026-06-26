/**
 * @file hybrid_examples.h
 * @brief Declarations for canonical hybrid automaton example constructors
 *
 * These functions construct standard hybrid automaton benchmarks
 * used in research and teaching (L6-L7).
 */

#ifndef HYBRID_EXAMPLES_H
#define HYBRID_EXAMPLES_H

#include "hybrid_automaton.h"
#include <stdbool.h>

HybridAutomaton* example_bouncing_ball(double h0, double v0, double c, double g);
HybridAutomaton* example_thermostat(double T_lo, double T_hi, double T_env,
                                     double K_h, double K_on,
                                     double T_init, bool start_on);
HybridAutomaton* example_two_tank(double K_inter, double K_out,
                                   double h1_init, double h2_init);
HybridAutomaton* example_train_gate(double approach_dist, double crossing_start,
                                     double crossing_end, double gate_down,
                                     double gate_up, double train_speed);
HybridAutomaton* example_dcdc_converter(double L, double C, double R,
                                         double V_in, double iL0, double vC0);
HybridAutomaton* example_engine_afr_control(void);
HybridAutomaton* example_robot_obstacle_avoidance(double cruise_speed,
                                                   double turn_radius,
                                                   double obstacle_x,
                                                   double obstacle_y);
HybridAutomaton* example_infusion_pump(void);

#endif /* HYBRID_EXAMPLES_H */
