CC=gcc
CFLAGS=-std=c11 -Wall -Wextra -O2 -Iinclude
LIB=libimpulsive.a
SRCS=$(wildcard src/*.c)
OBJS=$(SRCS:.c=.o)

all: $(LIB)

$(LIB): $(OBJS)
	ar rcs $@ $^

src/%.o: src/%.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@

test: $(LIB)
	for t in tests/test_*.c; do \
		$(CC) $(CFLAGS) $$t -L. -limpulsive -lm -o $$t.out && ./$$t.out; \
		rm -f $$t.out; \
	done

examples: $(LIB)
	for e in examples/example_*.c; do \
		$(CC) $(CFLAGS) $$e -L. -limpulsive -lm -o $$e.out && ./$$e.out; \
		rm -f $$e.out; \
	done

demos: $(LIB)
	for d in demos/demo_*.c; do \
		$(CC) $(CFLAGS) $$d -L. -limpulsive -lm -o $$d.out && ./$$d.out; \
		rm -f $$d.out; \
	done

bench: $(LIB)
	$(CC) $(CFLAGS) benches/bench_impulsive.c -L. -limpulsive -lm -o bench.out && ./bench.out; rm -f bench.out

clean:
	rm -f src/*.o $(LIB) tests/*.out examples/*.out demos/*.out bench.out

.PHONY: all test examples demos bench clean
