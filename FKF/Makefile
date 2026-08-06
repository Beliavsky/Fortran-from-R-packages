FC = gfortran
BUILD ?= build
MODE ?= checked

BASE_FLAGS = -std=f2018 -Wall -Wextra -Wpedantic -Werror -fimplicit-none -ffree-line-length-none
CHECKED_FLAGS = -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan
OPTIMIZED_FLAGS = -O3 -march=native
ifeq ($(MODE),optimized)
  FFLAGS = $(BASE_FLAGS) $(OPTIMIZED_FLAGS)
else
  FFLAGS = $(BASE_FLAGS) $(CHECKED_FLAGS)
endif

SRC = src/fkf_kinds.f90 src/fkf_linalg.f90 src/fkf_types.f90 src/fkf_filter.f90 src/fkf_smoother.f90 src/fkf.f90
OBJ = $(patsubst src/%.f90,$(BUILD)/%.o,$(SRC))
TESTS = test_scalar_filter test_missing_filter test_time_varying test_smoother test_likelihood_status

.PHONY: all test example clean checked optimized
all: test example

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/%.o: src/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/test_support.o: test/test_support.f90 $(OBJ)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/test_%: test/test_%.f90 $(OBJ) $(BUILD)/test_support.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $^ -o $@

$(BUILD)/local_level_example: example/local_level_example.f90 $(OBJ)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $^ -o $@

test: $(addprefix $(BUILD)/,$(TESTS))
	@set -e; for t in $(TESTS); do ./$(BUILD)/$$t; done

example: $(BUILD)/local_level_example
	./$(BUILD)/local_level_example

checked:
	$(MAKE) clean
	$(MAKE) MODE=checked all

optimized:
	$(MAKE) clean
	$(MAKE) MODE=optimized all

clean:
	rm -rf build
