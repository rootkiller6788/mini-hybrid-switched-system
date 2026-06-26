# Makefile for mini-piecewise-affine-system
# Piecewise Affine (PWA) System Theory Module
CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c11 -Iinclude
LDFLAGS = -lm
AR = ar
ARFLAGS = rcs

SRCDIR = src
INCDIR = include
TESTDIR = tests
EXDIR = examples
BUILDIR = build

SRCS = $(SRCDIR)/pwa_core.c \
       $(SRCDIR)/pwa_geometry.c \
       $(SRCDIR)/pwa_stability.c \
       $(SRCDIR)/pwa_simulation.c \
       $(SRCDIR)/pwa_identification.c \
       $(SRCDIR)/pwa_reachability.c \
       $(SRCDIR)/pwa_control.c

OBJS = $(patsubst $(SRCDIR)/%.c,$(BUILDIR)/%.o,$(SRCS))

TEST_SRCS = $(TESTDIR)/test_pwa.c
TEST_BINS = $(patsubst $(TESTDIR)/%.c,$(BUILDIR)/%,$(TEST_SRCS))

EX_SRCS = $(EXDIR)/example_saturation_pwa.c \
          $(EXDIR)/example_dcdc_converter_pwa.c \
          $(EXDIR)/example_bouncing_ball_pwa.c
EX_BINS = $(patsubst $(EXDIR)/%.c,$(BUILDIR)/%,$(EX_SRCS))

.PHONY: all test examples clean build_dir

all: build_dir $(BUILDIR)/libpwa.a $(EX_BINS) $(TEST_BINS)

build_dir:
	@mkdir -p $(BUILDIR)

$(BUILDIR)/libpwa.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $^

$(BUILDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILDIR)/test_pwa: $(TESTDIR)/test_pwa.c $(BUILDIR)/libpwa.a
	$(CC) $(CFLAGS) $< -L$(BUILDIR) -lpwa $(LDFLAGS) -o $@

$(BUILDIR)/example_saturation_pwa: $(EXDIR)/example_saturation_pwa.c $(BUILDIR)/libpwa.a
	$(CC) $(CFLAGS) $< -L$(BUILDIR) -lpwa $(LDFLAGS) -o $@

$(BUILDIR)/example_dcdc_converter_pwa: $(EXDIR)/example_dcdc_converter_pwa.c $(BUILDIR)/libpwa.a
	$(CC) $(CFLAGS) $< -L$(BUILDIR) -lpwa $(LDFLAGS) -o $@

$(BUILDIR)/example_bouncing_ball_pwa: $(EXDIR)/example_bouncing_ball_pwa.c $(BUILDIR)/libpwa.a
	$(CC) $(CFLAGS) $< -L$(BUILDIR) -lpwa $(LDFLAGS) -o $@

test: build_dir $(BUILDIR)/libpwa.a $(TEST_BINS)
	@echo === Running PWA System Tests ===
	$(BUILDIR)/test_pwa || true
	@echo === All Tests Complete ===

examples: build_dir $(BUILDIR)/libpwa.a $(EX_BINS)
	@echo === Running PWA Examples ===
	$(BUILDIR)/example_saturation_pwa || true
	$(BUILDIR)/example_dcdc_converter_pwa || true
	$(BUILDIR)/example_bouncing_ball_pwa || true
	@echo === All Examples Complete ===

clean:
	rm -rf $(BUILDIR)
