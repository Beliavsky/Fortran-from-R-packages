FC = gfortran
.NOTPARALLEL:
MODE ?= checked
BUILD_DIR := build/$(MODE)
MOD_DIR := $(BUILD_DIR)/mod
OBJ_DIR := $(BUILD_DIR)/obj
BIN_DIR := $(BUILD_DIR)/bin

BASE_FLAGS := -std=f2018 -ffree-line-length-none -J$(MOD_DIR) -I$(MOD_DIR)
CHECKED_FLAGS := -O0 -g -Wall -Wextra -Wpedantic -Wconversion-extra -fcheck=all -fbacktrace
OPTIMIZED_FLAGS := -O3 -march=native -Wall -Wextra -Wpedantic
ifeq ($(MODE),optimized)
  FFLAGS := $(BASE_FLAGS) $(OPTIMIZED_FLAGS)
else
  FFLAGS := $(BASE_FLAGS) $(CHECKED_FLAGS)
endif

SOURCES := \
  vendor/rpeif/vendor/robstattm-lite/src/robstattm_kinds.f90 \
  vendor/rpeif/vendor/robstattm-lite/src/robstattm_utils.f90 \
  vendor/rpeif/vendor/robstattm-lite/src/robstattm_psi.f90 \
  vendor/rpeif/src/rpeif_kinds.f90 \
  vendor/rpeif/src/rpeif_types.f90 \
  vendor/rpeif/src/rpeif_stats.f90 \
  vendor/rpeif/src/rpeif_robust.f90 \
  vendor/rpeif/src/rpeif_prewhiten.f90 \
  vendor/rpeif/src/rpeif_nuisance.f90 \
  vendor/rpeif/src/rpeif_influence.f90 \
  vendor/rpeif/src/rpeif_compat.f90 \
  vendor/rpeif/src/rpeif.f90 \
  vendor/rpeglmen/src/rpeglmen_kinds.f90 \
  vendor/rpeglmen/src/rpeglmen_types.f90 \
  vendor/rpeglmen/src/rpeglmen_math.f90 \
  vendor/rpeglmen/src/rpeglmen_penalty.f90 \
  vendor/rpeglmen/src/rpeglmen_likelihood.f90 \
  vendor/rpeglmen/src/rpeglmen_solver.f90 \
  vendor/rpeglmen/src/rpeglmen_gamma_mle.f90 \
  vendor/rpeglmen/src/rpeglmen_cv.f90 \
  vendor/rpeglmen/src/rpeglmen_compat.f90 \
  vendor/rpeglmen/src/rpeglmen.f90 \
  src/rpese_kinds.f90 \
  src/rpese_types.f90 \
  src/rpese_measures.f90 \
  src/rpese_timeseries.f90 \
  src/rpese_bootstrap.f90 \
  src/rpese_core.f90 \
  src/rpese_api.f90 \
  src/rpese_compat.f90 \
  src/rpese.f90

OBJECTS := $(patsubst %.f90,$(OBJ_DIR)/%.o,$(SOURCES))
TESTS := $(patsubst test/%.f90,$(BIN_DIR)/%,$(wildcard test/*.f90))
EXAMPLES := $(patsubst example/%.f90,$(BIN_DIR)/%,$(wildcard example/*.f90))
LIBRARY := $(BUILD_DIR)/librpese.a

.PHONY: all test examples clean
all: $(LIBRARY)

$(LIBRARY): $(OBJECTS)
	@mkdir -p $(dir $@)
	ar rcs $@ $(OBJECTS)

$(OBJ_DIR)/%.o: %.f90
	@mkdir -p $(dir $@) $(MOD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BIN_DIR)/%: test/%.f90 $(LIBRARY)
	@mkdir -p $(BIN_DIR)
	$(FC) $(FFLAGS) $< $(LIBRARY) -o $@

$(BIN_DIR)/%: example/%.f90 $(LIBRARY)
	@mkdir -p $(BIN_DIR)
	$(FC) $(FFLAGS) $< $(LIBRARY) -o $@

test: $(TESTS)
	@set -e; for program in $(TESTS); do echo "== $$program =="; $$program; done

examples: $(EXAMPLES)
	@set -e; for program in $(EXAMPLES); do echo "== $$program =="; $$program; done

clean:
	rm -rf build
