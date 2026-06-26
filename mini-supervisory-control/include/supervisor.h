/**
 * @file supervisor.h
 * @brief Supervisor Definition and Operations
 *
 * A supervisor S observes events and outputs control patterns.
 * Ramadge-Wonham framework: S: L(G) -> Gamma
 */

#ifndef SUPERVISOR_H
#define SUPERVISOR_H

#include "des_automaton.h"
#include <stddef.h>
#include <stdint.h>

#define SUP_MAX_PATTERNS   1024
#define SUP_MAX_HISTORY    512
#define SUP_MAX_MODULES    16

typedef struct {
    uint16_t state;
    uint32_t enabled_mask;
} sup_control_pattern_t;

typedef enum {
    SUP_TYPE_MONOLITHIC = 0,
    SUP_TYPE_MODULAR    = 1,
    SUP_TYPE_DECENTRALIZED = 2,
    SUP_TYPE_HIERARCHICAL = 3,
} sup_type_t;

/* Forward declaration to allow self-referential pointer */
typedef struct supervisor_t supervisor_t;

struct supervisor_t {
    sup_type_t          type;
    des_automaton_t    *plant;
    des_automaton_t    *spec;
    sup_control_pattern_t patterns[SUP_MAX_PATTERNS];
    uint16_t            npatterns;
    uint16_t            history[SUP_MAX_HISTORY];
    uint16_t            hist_len;
    supervisor_t       *modules[SUP_MAX_MODULES];
    uint16_t            nmodules;
    char                name[64];
    uint32_t            events_observed;
    uint32_t            events_disabled;
    uint32_t            events_enabled;
    uint32_t            violations;
    uint16_t            current_state;
};

void supervisor_init(supervisor_t *S, des_automaton_t *plant,
                      des_automaton_t *spec, sup_type_t type,
                      const char *name);
void supervisor_set_pattern(supervisor_t *S, uint16_t state,
                             uint32_t enabled);
uint32_t supervisor_current_pattern(const supervisor_t *S, uint16_t state);
int supervisor_is_enabled(const supervisor_t *S, uint16_t ev);
int supervisor_observe(supervisor_t *S, uint16_t ev);
uint16_t supervisor_enabled_set(const supervisor_t *S, uint16_t *enabled,
                                 uint16_t cap);
uint16_t supervisor_disabled_set(const supervisor_t *S, uint16_t *disabled,
                                  uint16_t cap);
int supervisor_is_proper(const supervisor_t *S);
void supervisor_closed_loop(const supervisor_t *S, des_automaton_t *closed_loop);
int supervisor_add_module(supervisor_t *S, supervisor_t *module_sup);
uint32_t supervisor_modular_pattern(const supervisor_t *S, uint16_t state);
void supervisor_reset(supervisor_t *S);
void supervisor_print_status(const supervisor_t *S, FILE *fp);
void supervisor_destroy(supervisor_t *S);

#endif /* SUPERVISOR_H */
