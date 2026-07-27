FC = gfortran
FFLAGS_COMMON = -std=f2018 -Wall -Wextra -Wimplicit-interface -Werror -ffree-line-length-none
FFLAGS_DEBUG = $(FFLAGS_COMMON) -O0 -g -fcheck=all -fbacktrace
FFLAGS_RELEASE = $(FFLAGS_COMMON) -O2 -fbacktrace
LIBS = -llapack -lblas
BUILD ?= build/debug
FFLAGS ?= $(FFLAGS_DEBUG)

SRC = src/pa_kinds.f90 src/pa_types.f90 src/pa_linalg.f90 src/pa_utils.f90 \
      src/pa_exposure.f90 src/pa_brinson.f90 src/pa_regression.f90 src/pa.f90
OBJ = $(patsubst src/%.f90,$(BUILD)/%.o,$(SRC))

.PHONY: all debug release test check apps clean

all: debug

debug:
	$(MAKE) BUILD=build/debug FFLAGS='$(FFLAGS_DEBUG)' check

release:
	$(MAKE) BUILD=build/release FFLAGS='$(FFLAGS_RELEASE)' check

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/%.o: src/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/test_pa: $(OBJ) test/test_pa.f90
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $(OBJ) test/test_pa.f90 $(LIBS) -o $@

$(BUILD)/demo_pa: $(OBJ) app/demo_pa.f90
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $(OBJ) app/demo_pa.f90 $(LIBS) -o $@

$(BUILD)/fit_csv: $(OBJ) app/fit_csv.f90
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $(OBJ) app/fit_csv.f90 $(LIBS) -o $@

$(BUILD)/attribution_example: $(OBJ) example/attribution_example.f90
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $(OBJ) example/attribution_example.f90 $(LIBS) -o $@

test: $(BUILD)/test_pa
	$(BUILD)/test_pa

apps: $(BUILD)/demo_pa $(BUILD)/fit_csv $(BUILD)/attribution_example
	$(BUILD)/demo_pa
	$(BUILD)/attribution_example
	$(BUILD)/fit_csv data/example_attribution.csv brinson geometric
	$(BUILD)/fit_csv data/example_attribution.csv regression geometric

check: test apps
	bash test/check_license.sh

clean:
	rm -rf build
