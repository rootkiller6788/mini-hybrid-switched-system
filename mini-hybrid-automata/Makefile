# Makefile for mini-hybrid-automata
#
# Builds the library, test suite, and example programs.
# Target: make test — compile and run all tests
#

CC      := gcc
CFLAGS  := -std=c11 -Wall -Wextra -O2 -g -Iinclude
LDFLAGS := -lm

BUILD_DIR := build
SRC_FILES := $(wildcard src/*.c)
OBJ_FILES := $(patsubst src/%.c,$(BUILD_DIR)/%.o,$(SRC_FILES))
TEST_SRC  := tests/test_hybrid_automata.c
TEST_BIN  := $(BUILD_DIR)/test_hybrid_automata
EXAMPLE_SRCS := $(wildcard examples/*.c)
EXAMPLE_BINS := $(patsubst examples/%.c,$(BUILD_DIR)/example_%,$(EXAMPLE_SRCS))

.PHONY: all test clean examples

all: $(TEST_BIN) examples

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: src/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Test binary links against all source object files
$(TEST_BIN): $(OBJ_FILES) $(TEST_SRC) | $(BUILD_DIR)
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) $(OBJ_FILES) $(LDFLAGS)

# Example binaries
$(BUILD_DIR)/example_%: examples/%.c $(OBJ_FILES) | $(BUILD_DIR)
	$(CC) $(CFLAGS) -o $@ $< $(OBJ_FILES) $(LDFLAGS)

# Build examples
examples: $(EXAMPLE_BINS)

# Build and run tests
test: $(TEST_BIN)
	@echo "=== Running tests ==="
	./$(TEST_BIN)

# Quick check: compile only
compile: $(OBJ_FILES)
	@echo "All source files compiled successfully."

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR)

# Print line counts
count:
	@echo "include/ headers: $$(wc -l include/*.h | tail -1)"
	@echo "src/ sources:     $$(wc -l src/*.c | tail -1)"
	@echo "Total include+src: $$(cat include/*.h src/*.c | wc -l)"
