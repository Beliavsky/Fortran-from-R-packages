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
  $(BUILD)/tvgarch_transition.o \
  $(BUILD)/tvgarch_model.o \
  $(BUILD)/tvgarch_multivariate.o \
  $(BUILD)/tvgarch_tests.o \
  $(BUILD)/tvgarch.o

TESTS = $(BUILD)/test_transition $(BUILD)/test_simulation $(BUILD)/test_fit \
        $(BUILD)/test_multivariate $(BUILD)/test_lm
APP = $(BUILD)/demo_tvgarch

.PHONY: all check optimized-check clean
all: $(TESTS) $(APP)

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
$(BUILD)/tvgarch_transition.o: src/tvgarch_transition.f90 $(BUILD)/garchx_kinds.o
	$(FC) $(FFLAGS) -J $(BUILD) -I $(BUILD) -c $< -o $@
$(BUILD)/tvgarch_model.o: src/tvgarch_model.f90 $(BUILD)/tvgarch_transition.o $(BUILD)/garchx_model.o
	$(FC) $(FFLAGS) -J $(BUILD) -I $(BUILD) -c $< -o $@
$(BUILD)/tvgarch_multivariate.o: src/tvgarch_multivariate.f90 $(BUILD)/tvgarch_model.o
	$(FC) $(FFLAGS) -J $(BUILD) -I $(BUILD) -c $< -o $@
$(BUILD)/tvgarch_tests.o: src/tvgarch_tests.f90 $(BUILD)/tvgarch_model.o
	$(FC) $(FFLAGS) -J $(BUILD) -I $(BUILD) -c $< -o $@
$(BUILD)/tvgarch.o: src/tvgarch.f90 $(BUILD)/tvgarch_multivariate.o $(BUILD)/tvgarch_tests.o $(BUILD)/garchx_inference.o
	$(FC) $(FFLAGS) -J $(BUILD) -I $(BUILD) -c $< -o $@

$(BUILD)/test_transition: test/test_transition.f90 $(MODS)
	$(FC) $(FFLAGS) -I $(BUILD) $< $(MODS) $(LDLIBS) -o $@
$(BUILD)/test_simulation: test/test_simulation.f90 $(MODS)
	$(FC) $(FFLAGS) -I $(BUILD) $< $(MODS) $(LDLIBS) -o $@
$(BUILD)/test_fit: test/test_fit.f90 $(MODS)
	$(FC) $(FFLAGS) -I $(BUILD) $< $(MODS) $(LDLIBS) -o $@
$(BUILD)/test_multivariate: test/test_multivariate.f90 $(MODS)
	$(FC) $(FFLAGS) -I $(BUILD) $< $(MODS) $(LDLIBS) -o $@
$(BUILD)/test_lm: test/test_lm.f90 $(MODS)
	$(FC) $(FFLAGS) -I $(BUILD) $< $(MODS) $(LDLIBS) -o $@
$(BUILD)/demo_tvgarch: app/demo_tvgarch.f90 $(MODS)
	$(FC) $(FFLAGS) -I $(BUILD) $< $(MODS) $(LDLIBS) -o $@

check: all
	@for t in $(TESTS); do $$t || exit 1; done
	$(APP) >/dev/null
	@echo "All tvgarch tests passed."

optimized-check:
	$(MAKE) clean
	$(MAKE) check BUILD=build/release FFLAGS='-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace'

clean:
	rm -rf build
