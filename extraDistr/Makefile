FC = gfortran
AR ?= ar
MODE ?= check
BUILD_DIR := build/$(MODE)
MOD_DIR := $(BUILD_DIR)/mod

BASE_FLAGS := -std=f2018 -ffree-line-length-none -fimplicit-none -Wall -Wextra -Werror \
              -Wno-compare-reals -Wno-integer-division
CHECK_FLAGS := -O0 -g -fcheck=all -fbacktrace
RELEASE_FLAGS := -O3

ifeq ($(MODE),release)
  FFLAGS := $(BASE_FLAGS) $(RELEASE_FLAGS)
else
  FFLAGS := $(BASE_FLAGS) $(CHECK_FLAGS)
endif

SRC_DIR := src
TEST_DIR := test
EXAMPLE_DIR := example

MODULES := extra_distr_kinds extra_distr_math extra_distr_rng \
           extra_distr_continuous extra_distr_discrete \
           extra_distr_multivariate extra_distr
OBJECTS := $(addprefix $(BUILD_DIR)/,$(addsuffix .o,$(MODULES)))
LIBRARY := $(BUILD_DIR)/libextra_distr.a
TEST_SUPPORT := $(BUILD_DIR)/test_support.o
TEST_NAMES := test_continuous test_discrete test_multivariate test_random test_api
TEST_EXES := $(addprefix $(BUILD_DIR)/,$(TEST_NAMES))
EXAMPLE_EXE := $(BUILD_DIR)/demo_extra_distr

.PHONY: all library check release run-tests example clean distclean

all: library
library: $(LIBRARY)

check:
	$(MAKE) MODE=check run-tests

release:
	$(MAKE) MODE=release run-tests

run-tests: $(TEST_EXES)
	@set -e; for t in $(TEST_EXES); do echo "Running $$t"; "$$t"; done

example: $(EXAMPLE_EXE)
	$(EXAMPLE_EXE)

$(BUILD_DIR) $(MOD_DIR):
	mkdir -p $@

$(BUILD_DIR)/extra_distr_kinds.o: $(SRC_DIR)/extra_distr_kinds.f90 | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(BUILD_DIR)/extra_distr_math.o: $(SRC_DIR)/extra_distr_math.f90 $(BUILD_DIR)/extra_distr_kinds.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(BUILD_DIR)/extra_distr_rng.o: $(SRC_DIR)/extra_distr_rng.f90 $(BUILD_DIR)/extra_distr_kinds.o $(BUILD_DIR)/extra_distr_math.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(BUILD_DIR)/extra_distr_continuous.o: $(SRC_DIR)/extra_distr_continuous.f90 $(BUILD_DIR)/extra_distr_kinds.o $(BUILD_DIR)/extra_distr_math.o $(BUILD_DIR)/extra_distr_rng.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(BUILD_DIR)/extra_distr_discrete.o: $(SRC_DIR)/extra_distr_discrete.f90 $(BUILD_DIR)/extra_distr_kinds.o $(BUILD_DIR)/extra_distr_math.o $(BUILD_DIR)/extra_distr_rng.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(BUILD_DIR)/extra_distr_multivariate.o: $(SRC_DIR)/extra_distr_multivariate.f90 $(BUILD_DIR)/extra_distr_kinds.o $(BUILD_DIR)/extra_distr_math.o $(BUILD_DIR)/extra_distr_rng.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(BUILD_DIR)/extra_distr.o: $(SRC_DIR)/extra_distr.f90 $(BUILD_DIR)/extra_distr_continuous.o $(BUILD_DIR)/extra_distr_discrete.o $(BUILD_DIR)/extra_distr_multivariate.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(LIBRARY): $(OBJECTS)
	$(AR) rcs $@ $^

$(TEST_SUPPORT): $(TEST_DIR)/test_support.f90 $(BUILD_DIR)/extra_distr_kinds.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(BUILD_DIR)/test_%: $(TEST_DIR)/test_%.f90 $(TEST_SUPPORT) $(LIBRARY) | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $< $(TEST_SUPPORT) $(LIBRARY) -o $@

$(EXAMPLE_EXE): $(EXAMPLE_DIR)/demo_extra_distr.f90 $(LIBRARY) | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $< $(LIBRARY) -o $@

clean:
	rm -rf build/check build/release

distclean: clean
	rm -f *.o *.mod
