FC = gfortran
BUILD ?= build/debug
FFLAGS ?= -std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace
LDLIBS ?= -llapack -lblas

MODS = \
  $(BUILD)/garchx_kinds.o \
  $(BUILD)/garchx_linalg.o \
  $(BUILD)/garchx_math.o \
  $(BUILD)/garchx_optimize.o \
  $(BUILD)/garchx_model.o \
  $(BUILD)/garchx_utilities.o \
  $(BUILD)/garchx_inference.o \
  $(BUILD)/garchx_csv.o

TESTS = $(BUILD)/test_core $(BUILD)/test_fit $(BUILD)/test_inference
APPS = $(BUILD)/demo_garchx $(BUILD)/fit_csv $(BUILD)/sparse_garchx_example

.PHONY: all check clean optimized-check license-check
all: $(TESTS) $(APPS)

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/garchx_kinds.o: src/garchx_kinds.f90 | $(BUILD)
	$(FC) $(FFLAGS) -J $(BUILD) -I $(BUILD) -c $< -o $@
$(BUILD)/garchx_linalg.o: src/garchx_linalg.f90 $(BUILD)/garchx_kinds.o
	$(FC) $(FFLAGS) -J $(BUILD) -I $(BUILD) -c $< -o $@
$(BUILD)/garchx_math.o: src/garchx_math.f90 $(BUILD)/garchx_kinds.o $(BUILD)/garchx_linalg.o
	$(FC) $(FFLAGS) -J $(BUILD) -I $(BUILD) -c $< -o $@
$(BUILD)/garchx_optimize.o: src/garchx_optimize.f90 $(BUILD)/garchx_kinds.o
	$(FC) $(FFLAGS) -J $(BUILD) -I $(BUILD) -c $< -o $@
$(BUILD)/garchx_model.o: src/garchx_model.f90 $(BUILD)/garchx_math.o $(BUILD)/garchx_optimize.o
	$(FC) $(FFLAGS) -J $(BUILD) -I $(BUILD) -c $< -o $@
$(BUILD)/garchx_utilities.o: src/garchx_utilities.f90 $(BUILD)/garchx_kinds.o
	$(FC) $(FFLAGS) -J $(BUILD) -I $(BUILD) -c $< -o $@
$(BUILD)/garchx_inference.o: src/garchx_inference.f90 $(BUILD)/garchx_math.o
	$(FC) $(FFLAGS) -J $(BUILD) -I $(BUILD) -c $< -o $@
$(BUILD)/garchx_csv.o: src/garchx_csv.f90 $(BUILD)/garchx_kinds.o
	$(FC) $(FFLAGS) -J $(BUILD) -I $(BUILD) -c $< -o $@

$(BUILD)/test_core: test/test_core.f90 $(MODS)
	$(FC) $(FFLAGS) -I $(BUILD) $< $(MODS) $(LDLIBS) -o $@
$(BUILD)/test_fit: test/test_fit.f90 $(MODS)
	$(FC) $(FFLAGS) -I $(BUILD) $< $(MODS) $(LDLIBS) -o $@
$(BUILD)/test_inference: test/test_inference.f90 $(MODS)
	$(FC) $(FFLAGS) -I $(BUILD) $< $(MODS) $(LDLIBS) -o $@
$(BUILD)/demo_garchx: app/demo_garchx.f90 $(MODS)
	$(FC) $(FFLAGS) -I $(BUILD) $< $(MODS) $(LDLIBS) -o $@
$(BUILD)/fit_csv: app/fit_csv.f90 $(MODS)
	$(FC) $(FFLAGS) -I $(BUILD) $< $(MODS) $(LDLIBS) -o $@
$(BUILD)/sparse_garchx_example: example/sparse_garchx_example.f90 $(MODS)
	$(FC) $(FFLAGS) -I $(BUILD) $< $(MODS) $(LDLIBS) -o $@

license-check:
	./scripts/check_license.sh

check: license-check all
	$(BUILD)/test_core
	$(BUILD)/test_fit
	$(BUILD)/test_inference
	$(BUILD)/demo_garchx >/dev/null
	$(BUILD)/sparse_garchx_example >/dev/null
	$(BUILD)/fit_csv data/example.csv 1 1 1 ordinary >/dev/null
	$(BUILD)/fit_csv data/example.csv 1 1 1 robust >/dev/null
	$(BUILD)/fit_csv data/example.csv 1 1 1 hac >/dev/null
	@echo "Application and CSV workflow tests passed."

optimized-check:
	$(MAKE) clean
	$(MAKE) check BUILD=build/release FFLAGS='-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace'

clean:
	rm -rf build
