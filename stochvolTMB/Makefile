FC = gfortran
AR ?= ar
MODE ?= checked

SRC = \
  src/stochvoltmb_kinds.f90 \
  src/stochvoltmb_status.f90 \
  src/stochvoltmb_rng.f90 \
  src/stochvoltmb_math.f90 \
  src/stochvoltmb_linalg.f90 \
  src/stochvoltmb_types.f90 \
  src/stochvoltmb_model.f90 \
  src/stochvoltmb_simulation.f90 \
  src/stochvoltmb.f90

TEST_SRC = $(wildcard test/test_*.f90)
TEST_BIN = $(patsubst test/%.f90,build/bin/%,$(TEST_SRC))
OBJ = $(patsubst src/%.f90,build/obj/%.o,$(SRC))

BASE_FLAGS = -std=f2018 -Wall -Wextra -Werror -Jbuild/mod -Ibuild/mod
ifeq ($(MODE),optimized)
  FFLAGS = $(BASE_FLAGS) -O3 -march=native
else
  FFLAGS = $(BASE_FLAGS) -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
endif

.PHONY: all test demo clean checked optimized
all: build/libstochvoltmb.a
checked:
	$(MAKE) clean
	$(MAKE) MODE=checked test
optimized:
	$(MAKE) clean
	$(MAKE) MODE=optimized test

test: $(TEST_BIN)
	@set -e; for t in $(TEST_BIN); do echo "Running $$t"; $$t; done

demo: build/bin/demo_stochvoltmb
	build/bin/demo_stochvoltmb

build/libstochvoltmb.a: $(OBJ)
	$(AR) rcs $@ $(OBJ)

build/obj/%.o: src/%.f90
	@mkdir -p build/obj build/mod build/bin
	$(FC) $(FFLAGS) -c $< -o $@

build/bin/test_%: test/test_%.f90 build/libstochvoltmb.a
	@mkdir -p build/bin
	$(FC) $(FFLAGS) $< build/libstochvoltmb.a -o $@

build/bin/demo_stochvoltmb: app/demo_stochvoltmb.f90 build/libstochvoltmb.a
	@mkdir -p build/bin
	$(FC) $(FFLAGS) $< build/libstochvoltmb.a -o $@

clean:
	rm -rf build *.mod *.o
