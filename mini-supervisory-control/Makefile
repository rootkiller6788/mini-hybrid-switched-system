# Makefile for mini-supervisory-control
CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS  = -lm

SRC = src/des_automaton.c src/supervisor.c src/controllability.c src/synthesis.c src/hybrid_supervisor.c
TEST_SRC = tests/test_supervisory.c
TEST_BIN = build/test_supervisory

EXAMPLE_SRC = examples/example_transfer_line.c examples/example_small_factory.c examples/example_switched_system.c

.PHONY: all test examples clean

all: build_dir test examples

build_dir:
	mkdir -p build

test: build_dir $(TEST_BIN)
	./$(TEST_BIN)

$(TEST_BIN): $(SRC) $(TEST_SRC)
	$(CC) $(CFLAGS) -o $@ $(SRC) $(TEST_SRC) $(LDFLAGS)

examples: build_dir build/example_transfer_line build/example_small_factory build/example_switched_system
	./build/example_transfer_line
	./build/example_small_factory
	./build/example_switched_system

build/example_transfer_line: $(SRC) examples/example_transfer_line.c
	$(CC) $(CFLAGS) -o $@ $(SRC) examples/example_transfer_line.c $(LDFLAGS)

build/example_small_factory: $(SRC) examples/example_small_factory.c
	$(CC) $(CFLAGS) -o $@ $(SRC) examples/example_small_factory.c $(LDFLAGS)

build/example_switched_system: $(SRC) examples/example_switched_system.c
	$(CC) $(CFLAGS) -o $@ $(SRC) examples/example_switched_system.c $(LDFLAGS)

lines:
	@echo include/ + src/ lines:
	@cat include/*.h src/*.c | wc -l

clean:
	rm -rf build/
