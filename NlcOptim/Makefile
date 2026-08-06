FC = gfortran
BUILD ?= build/check
MODE ?= check

BASE_FLAGS = -std=f2018 -Wall -Wextra -Werror -pedantic -Wno-compare-reals
CHECK_FLAGS = -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
OPT_FLAGS = -O3
FFLAGS = $(BASE_FLAGS) $(if $(filter $(MODE),check),$(CHECK_FLAGS),$(OPT_FLAGS))

SRC = src/quadprog_kinds.f90 src/quadprog_core.f90 src/quadprog.f90 src/nlcoptim.f90
OBJ = $(BUILD)/quadprog_kinds.o $(BUILD)/quadprog_core.o $(BUILD)/quadprog.o $(BUILD)/nlcoptim.o
TESTS = test_unconstrained test_linear_constraints test_nonlinear_equality test_nonlinear_inequality test_qp_and_errors

.PHONY: all check optimized test example clean

all: check

check:
	$(MAKE) MODE=check BUILD=build/check test example

optimized:
	$(MAKE) MODE=optimized BUILD=build/optimized test example

test: $(addprefix $(BUILD)/,$(TESTS))
	@set -e; for t in $(TESTS); do ./$(BUILD)/$$t; done

example: $(BUILD)/demo_nlcoptim
	./$(BUILD)/demo_nlcoptim

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/quadprog_kinds.o: src/quadprog_kinds.f90 | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/quadprog_core.o: src/quadprog_core.f90 $(BUILD)/quadprog_kinds.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/quadprog.o: src/quadprog.f90 $(BUILD)/quadprog_core.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/nlcoptim.o: src/nlcoptim.f90 $(BUILD)/quadprog.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/%: test/%.f90 $(OBJ)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $< $(OBJ) -o $@

$(BUILD)/demo_nlcoptim: app/demo_nlcoptim.f90 $(OBJ)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $< $(OBJ) -o $@

clean:
	rm -rf build
