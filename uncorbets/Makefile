FC = gfortran
BUILD ?= build/check
MODE ?= check

BASE_FLAGS = -std=f2018 -Wall -Wextra -Werror -pedantic
CHECK_FLAGS = -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
OPT_FLAGS = -O3
FFLAGS = $(BASE_FLAGS) $(if $(filter $(MODE),check),$(CHECK_FLAGS),$(OPT_FLAGS))

SRC = src/uncorbets_kinds.f90 src/uncorbets_types.f90 src/uncorbets_linalg.f90 \
      src/uncorbets_core.f90 src/uncorbets.f90
OBJ = $(BUILD)/uncorbets_kinds.o $(BUILD)/uncorbets_types.o \
      $(BUILD)/uncorbets_linalg.o $(BUILD)/uncorbets_core.o $(BUILD)/uncorbets.o
TESTS = test_sqrtm test_torsion test_effective_bets test_max_effective_bets test_errors

.PHONY: all check optimized test app example clean

all: check

check:
	$(MAKE) MODE=check BUILD=build/check test app example

optimized:
	$(MAKE) MODE=optimized BUILD=build/optimized test app example

test: $(addprefix $(BUILD)/,$(TESTS))
	@set -e; for t in $(TESTS); do ./$(BUILD)/$$t; done

app: $(BUILD)/demo_uncorbets
	./$(BUILD)/demo_uncorbets

example: $(BUILD)/basic_uncorbets
	./$(BUILD)/basic_uncorbets

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/uncorbets_kinds.o: src/uncorbets_kinds.f90 | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/uncorbets_types.o: src/uncorbets_types.f90 $(BUILD)/uncorbets_kinds.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/uncorbets_linalg.o: src/uncorbets_linalg.f90 $(BUILD)/uncorbets_types.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/uncorbets_core.o: src/uncorbets_core.f90 $(BUILD)/uncorbets_linalg.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/uncorbets.o: src/uncorbets.f90 $(BUILD)/uncorbets_core.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/%: test/%.f90 $(OBJ)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $< $(OBJ) -o $@

$(BUILD)/demo_uncorbets: app/demo_uncorbets.f90 $(OBJ)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $< $(OBJ) -o $@

$(BUILD)/basic_uncorbets: example/basic_uncorbets.f90 $(OBJ)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $< $(OBJ) -o $@

clean:
	rm -rf build
