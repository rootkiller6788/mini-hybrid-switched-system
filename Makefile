# Makefile for mini-hybrid-switched-system (Module 19)
# Top-level integration: Hybrid Switched System Unified Framework
#
# This Makefile builds the unified HSS library and its tests.
# The sub-modules (mini-dwell-time-analysis, mini-hybrid-automata, etc.)
# each have their own Makefiles and are built independently.

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -pedantic -O2 -g
CFLAGS  += -Iinclude

# HSS integration library
LIB_SRC  = src/hss_core.c \
           src/hss_simulation.c \
           src/hss_analysis.c \
           src/hss_hybrid_stability.c \
           src/hss_integrator.c

LIB_OBJ  = build/hss_core.o \
           build/hss_simulation.o \
           build/hss_analysis.o \
           build/hss_hybrid_stability.o \
           build/hss_integrator.o

LIB      = build/libhss.a

# Test sources
TEST_SRC = tests/test_hss_core.c \
           tests/test_hss_simulation.c \
           tests/test_hss_analysis.c \
           tests/test_hss_stability.c \
           tests/test_hss_integrator.c

TEST_BIN = build/test_hss_core \
           build/test_hss_simulation \
           build/test_hss_analysis \
           build/test_hss_stability \
           build/test_hss_integrator

# Example binaries
EXAMPLES = examples/example_bouncing_ball \
           examples/example_thermostat \
           examples/example_dcdc_converter \
           examples/example_cruise_control

# Benchmarks
BENCHES  = benches/bench_hss_sim

.PHONY: all clean test lib build-dir examples benches run-tests

all: lib test examples

build-dir:
	@mkdir -p build

# ---- Library ----
lib: build-dir $(LIB)

$(LIB): $(LIB_OBJ)
	ar rcs $(LIB) $(LIB_OBJ)

build/hss_core.o: src/hss_core.c include/hss_core.h
	$(CC) $(CFLAGS) -c src/hss_core.c -o build/hss_core.o

build/hss_simulation.o: src/hss_simulation.c include/hss_simulation.h include/hss_core.h
	$(CC) $(CFLAGS) -c src/hss_simulation.c -o build/hss_simulation.o

build/hss_analysis.o: src/hss_analysis.c include/hss_analysis.h include/hss_core.h include/hss_simulation.h
	$(CC) $(CFLAGS) -c src/hss_analysis.c -o build/hss_analysis.o

build/hss_hybrid_stability.o: src/hss_hybrid_stability.c include/hss_hybrid_stability.h include/hss_core.h include/hss_analysis.h
	$(CC) $(CFLAGS) -c src/hss_hybrid_stability.c -o build/hss_hybrid_stability.o

build/hss_integrator.o: src/hss_integrator.c include/hss_integrator.h include/hss_core.h include/hss_analysis.h
	$(CC) $(CFLAGS) -c src/hss_integrator.c -o build/hss_integrator.o

# ---- Tests ----
test: build-dir lib $(TEST_BIN)
	@echo "============================================"
	@echo "Running Hybrid Switched System Tests"
	@echo "============================================"
	@FAILED=0; \
	for t in $(TEST_BIN); do \
		echo "--- Running $$t ---"; \
		./$$t && echo "PASS: $$t" || { echo "FAIL: $$t"; FAILED=1; }; \
		echo ""; \
	done; \
	if [ $$FAILED -eq 0 ]; then \
		echo "=== All HSS tests passed ==="; \
	else \
		echo "=== Some HSS tests FAILED ==="; \
		exit 1; \
	fi

build/test_hss_core: tests/test_hss_core.c $(LIB)
	$(CC) $(CFLAGS) -o build/test_hss_core tests/test_hss_core.c -Lbuild -lhss -lm

build/test_hss_simulation: tests/test_hss_simulation.c $(LIB)
	$(CC) $(CFLAGS) -o build/test_hss_simulation tests/test_hss_simulation.c -Lbuild -lhss -lm

build/test_hss_analysis: tests/test_hss_analysis.c $(LIB)
	$(CC) $(CFLAGS) -o build/test_hss_analysis tests/test_hss_analysis.c -Lbuild -lhss -lm

build/test_hss_stability: tests/test_hss_stability.c $(LIB)
	$(CC) $(CFLAGS) -o build/test_hss_stability tests/test_hss_stability.c -Lbuild -lhss -lm

build/test_hss_integrator: tests/test_hss_integrator.c $(LIB)
	$(CC) $(CFLAGS) -o build/test_hss_integrator tests/test_hss_integrator.c -Lbuild -lhss -lm

# ---- Examples ----
examples: build-dir lib $(EXAMPLES)

examples/example_bouncing_ball: examples/example_bouncing_ball.c $(LIB)
	$(CC) $(CFLAGS) -o examples/example_bouncing_ball examples/example_bouncing_ball.c -Lbuild -lhss -lm

examples/example_thermostat: examples/example_thermostat.c $(LIB)
	$(CC) $(CFLAGS) -o examples/example_thermostat examples/example_thermostat.c -Lbuild -lhss -lm

examples/example_dcdc_converter: examples/example_dcdc_converter.c $(LIB)
	$(CC) $(CFLAGS) -o examples/example_dcdc_converter examples/example_dcdc_converter.c -Lbuild -lhss -lm

examples/example_cruise_control: examples/example_cruise_control.c $(LIB)
	$(CC) $(CFLAGS) -o examples/example_cruise_control examples/example_cruise_control.c -Lbuild -lhss -lm

# ---- Benchmarks ----
benches: build-dir lib $(BENCHES)

benches/bench_hss_sim: benches/bench_hss_sim.c $(LIB)
	$(CC) $(CFLAGS) -o benches/bench_hss_sim benches/bench_hss_sim.c -Lbuild -lhss -lm

# ---- All sub-modules test ----
test-all: test
	@echo ""
	@echo "=== Testing sub-modules ==="
	@for d in mini-*/; do \
		echo "--- $$d ---"; \
		$(MAKE) -C $$d test 2>&1 | tail -3 || true; \
		echo ""; \
	done

# ---- Clean ----
clean:
	rm -rf build
	rm -f examples/example_bouncing_ball
	rm -f examples/example_thermostat
	rm -f examples/example_dcdc_converter
	rm -f examples/example_cruise_control
	rm -f benches/bench_hss_sim

# ---- Help ----
help:
	@echo "Targets:"
	@echo "  all       - Build library, tests, and examples"
	@echo "  lib       - Build libhss.a"
	@echo "  test      - Build and run all tests"
	@echo "  test-all  - Run HSS tests + all sub-module tests"
	@echo "  examples  - Build example programs"
	@echo "  benches   - Build benchmarks"
	@echo "  clean     - Remove build artifacts"
