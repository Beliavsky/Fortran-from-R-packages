FC = gfortran
.NOTPARALLEL:
MODE ?= checked
BUILD_DIR := build/$(MODE)
MOD_DIR := $(BUILD_DIR)/mod
OBJ_DIR := $(BUILD_DIR)/obj
BIN_DIR := $(BUILD_DIR)/bin

BASE_FLAGS := -std=f2018 -ffree-line-length-none -J$(MOD_DIR) -I$(MOD_DIR)
ifeq ($(MODE),checked)
  FFLAGS := $(BASE_FLAGS) -O0 -g -Wall -Wextra -Wpedantic -Wconversion-extra -fcheck=all -fbacktrace
else
  FFLAGS := $(BASE_FLAGS) -O3 -Wall -Wextra -Wpedantic
endif

SOURCES := \
  vendor/robstattm-lite/src/robstattm_kinds.f90 \
  vendor/robstattm-lite/src/robstattm_utils.f90 \
  vendor/robstattm-lite/src/robstattm_psi.f90 \
  src/rpeif_kinds.f90 \
  src/rpeif_types.f90 \
  src/rpeif_stats.f90 \
  src/rpeif_robust.f90 \
  src/rpeif_prewhiten.f90 \
  src/rpeif_nuisance.f90 \
  src/rpeif_influence.f90 \
  src/rpeif_compat.f90 \
  src/rpeif.f90

OBJECTS := $(patsubst %.f90,$(OBJ_DIR)/%.o,$(SOURCES))
TESTS := $(patsubst test/%.f90,$(BIN_DIR)/%,$(wildcard test/*.f90))
EXAMPLES := $(patsubst example/%.f90,$(BIN_DIR)/%,$(wildcard example/*.f90))
LIBRARY := $(BUILD_DIR)/librpeif.a

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

clean:
	rm -rf build
