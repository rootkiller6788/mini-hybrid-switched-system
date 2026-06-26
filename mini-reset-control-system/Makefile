# Makefile for mini-reset-control-system

CC      ?= gcc
CFLAGS  ?= -std=c11 -Wall -Wextra -O2 -g -D_USE_MATH_DEFINES -DM_PI=3.14159265358979323846
AR      ?= ar
ARFLAGS ?= rcs

INC_DIR  = include
SRC_DIR  = src
BLD_DIR  = build
INCLUDES = -I$(INC_DIR)
LDFLAGS  = -lm

SRCS = src/reset_core.c src/reset_element.c src/reset_system.c src/reset_simulation.c src/reset_math.c src/reset_lyapunov.c
OBJS = $(SRCS:src/%.c=$(BLD_DIR)/%.o)
LIB  = $(BLD_DIR)/libreset.a

.PHONY: all test clean examples

all: $(LIB)

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

$(LIB): $(OBJS) | $(BLD_DIR)
	$(AR) $(ARFLAGS) $@ $^

$(BLD_DIR)/%.o: src/%.c | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

test: $(LIB)
	$(CC) $(CFLAGS) $(INCLUDES) tests/test_reset.c -L$(BLD_DIR) -lreset $(LDFLAGS) -o $(BLD_DIR)/test_reset
	./$(BLD_DIR)/test_reset
	$(CC) $(CFLAGS) $(INCLUDES) tests/test_elements.c -L$(BLD_DIR) -lreset $(LDFLAGS) -o $(BLD_DIR)/test_elements
	./$(BLD_DIR)/test_elements

examples: $(LIB)
	$(CC) $(CFLAGS) $(INCLUDES) examples/ex_clegg_step.c -L$(BLD_DIR) -lreset $(LDFLAGS) -o $(BLD_DIR)/ex_clegg_step
	./$(BLD_DIR)/ex_clegg_step
	$(CC) $(CFLAGS) $(INCLUDES) examples/ex_fore_response.c -L$(BLD_DIR) -lreset $(LDFLAGS) -o $(BLD_DIR)/ex_fore_response
	./$(BLD_DIR)/ex_fore_response
	$(CC) $(CFLAGS) $(INCLUDES) examples/ex_reset_pid.c -L$(BLD_DIR) -lreset $(LDFLAGS) -o $(BLD_DIR)/ex_reset_pid
	./$(BLD_DIR)/ex_reset_pid
	$(CC) $(CFLAGS) $(INCLUDES) examples/ex_reset_feedback.c -L$(BLD_DIR) -lreset $(LDFLAGS) -o $(BLD_DIR)/ex_reset_feedback
	./$(BLD_DIR)/ex_reset_feedback

clean:
	rm -rf $(BLD_DIR)
