FC = gfortran

BUILD_DIR = build
MOD_DIR = $(BUILD_DIR)/mod
OBJ_DIR = $(BUILD_DIR)/obj
BIN_DIR = $(BUILD_DIR)/bin

COMMON_WARNINGS = -std=f2018 -Wall -Wextra -Wimplicit-interface -Werror
FFLAGS ?= -O2 $(COMMON_WARNINGS)
MODULE_FLAGS = -J$(MOD_DIR) -I$(MOD_DIR)

LIB_OBJECTS = \
	$(OBJ_DIR)/bayesgarch_kinds.o \
	$(OBJ_DIR)/bayesgarch_rng.o \
	$(OBJ_DIR)/bayesgarch_math.o \
	$(OBJ_DIR)/bayesgarch_garch.o \
	$(OBJ_DIR)/bayesgarch_sampler.o \
	$(OBJ_DIR)/bayesgarch_sample.o

TEST_SUPPORT = $(OBJ_DIR)/test_support.o

.PHONY: all apps tests examples check optimized-check clean

all: apps

apps: $(BIN_DIR)/demo_bayesgarch $(BIN_DIR)/fit_csv

tests: $(BIN_DIR)/test_filters $(BIN_DIR)/test_rng_math $(BIN_DIR)/test_sampler

examples: $(BIN_DIR)/custom_constraint_example

$(MOD_DIR) $(OBJ_DIR) $(BIN_DIR):
	mkdir -p $@

$(OBJ_DIR)/bayesgarch_kinds.o: src/bayesgarch_kinds.f90 | $(MOD_DIR) $(OBJ_DIR)
	$(FC) $(FFLAGS) $(MODULE_FLAGS) -c $< -o $@

$(OBJ_DIR)/bayesgarch_rng.o: src/bayesgarch_rng.f90 $(OBJ_DIR)/bayesgarch_kinds.o | $(MOD_DIR) $(OBJ_DIR)
	$(FC) $(FFLAGS) $(MODULE_FLAGS) -c $< -o $@

$(OBJ_DIR)/bayesgarch_math.o: src/bayesgarch_math.f90 $(OBJ_DIR)/bayesgarch_kinds.o $(OBJ_DIR)/bayesgarch_rng.o | $(MOD_DIR) $(OBJ_DIR)
	$(FC) $(FFLAGS) $(MODULE_FLAGS) -c $< -o $@

$(OBJ_DIR)/bayesgarch_garch.o: src/bayesgarch_garch.f90 $(OBJ_DIR)/bayesgarch_kinds.o \
	$(OBJ_DIR)/bayesgarch_rng.o $(OBJ_DIR)/bayesgarch_math.o | $(MOD_DIR) $(OBJ_DIR)
	$(FC) $(FFLAGS) $(MODULE_FLAGS) -c $< -o $@

$(OBJ_DIR)/bayesgarch_sampler.o: src/bayesgarch_sampler.f90 $(OBJ_DIR)/bayesgarch_kinds.o \
	$(OBJ_DIR)/bayesgarch_rng.o $(OBJ_DIR)/bayesgarch_math.o $(OBJ_DIR)/bayesgarch_garch.o | $(MOD_DIR) $(OBJ_DIR)
	$(FC) $(FFLAGS) $(MODULE_FLAGS) -c $< -o $@

$(OBJ_DIR)/bayesgarch_sample.o: src/bayesgarch_sample.f90 $(OBJ_DIR)/bayesgarch_kinds.o \
	$(OBJ_DIR)/bayesgarch_sampler.o | $(MOD_DIR) $(OBJ_DIR)
	$(FC) $(FFLAGS) $(MODULE_FLAGS) -c $< -o $@

$(TEST_SUPPORT): test/test_support.f90 $(OBJ_DIR)/bayesgarch_kinds.o | $(MOD_DIR) $(OBJ_DIR)
	$(FC) $(FFLAGS) $(MODULE_FLAGS) -c $< -o $@

$(BIN_DIR)/demo_bayesgarch: app/demo_bayesgarch.f90 $(LIB_OBJECTS) | $(BIN_DIR)
	$(FC) $(FFLAGS) $(MODULE_FLAGS) $< $(LIB_OBJECTS) -o $@

$(BIN_DIR)/fit_csv: app/fit_csv.f90 $(LIB_OBJECTS) | $(BIN_DIR)
	$(FC) $(FFLAGS) $(MODULE_FLAGS) $< $(LIB_OBJECTS) -o $@

$(BIN_DIR)/custom_constraint_example: example/custom_constraint.f90 $(LIB_OBJECTS) | $(BIN_DIR)
	$(FC) $(FFLAGS) $(MODULE_FLAGS) $< $(LIB_OBJECTS) -o $@

$(BIN_DIR)/test_filters: test/test_filters.f90 $(LIB_OBJECTS) $(TEST_SUPPORT) | $(BIN_DIR)
	$(FC) $(FFLAGS) $(MODULE_FLAGS) $< $(LIB_OBJECTS) $(TEST_SUPPORT) -o $@

$(BIN_DIR)/test_rng_math: test/test_rng_math.f90 $(LIB_OBJECTS) $(TEST_SUPPORT) | $(BIN_DIR)
	$(FC) $(FFLAGS) $(MODULE_FLAGS) $< $(LIB_OBJECTS) $(TEST_SUPPORT) -o $@

$(BIN_DIR)/test_sampler: test/test_sampler.f90 $(LIB_OBJECTS) $(TEST_SUPPORT) | $(BIN_DIR)
	$(FC) $(FFLAGS) $(MODULE_FLAGS) $< $(LIB_OBJECTS) $(TEST_SUPPORT) -o $@

check:
	$(MAKE) clean
	$(MAKE) FFLAGS="-O0 -g $(COMMON_WARNINGS) -fcheck=all -fbacktrace" apps tests examples
	./scripts/check_license.sh
	$(BIN_DIR)/test_filters
	$(BIN_DIR)/test_rng_math
	$(BIN_DIR)/test_sampler
	$(BIN_DIR)/demo_bayesgarch > $(BUILD_DIR)/demo_output.txt
	$(BIN_DIR)/fit_csv data/example_returns.csv 250 50 5 1 > $(BUILD_DIR)/fit_csv_output.txt
	$(BIN_DIR)/custom_constraint_example > $(BUILD_DIR)/constraint_output.txt

optimized-check:
	$(MAKE) clean
	$(MAKE) FFLAGS="-O2 $(COMMON_WARNINGS)" apps tests examples
	./scripts/check_license.sh
	$(BIN_DIR)/test_filters
	$(BIN_DIR)/test_rng_math
	$(BIN_DIR)/test_sampler
	$(BIN_DIR)/demo_bayesgarch > $(BUILD_DIR)/demo_output.txt
	$(BIN_DIR)/fit_csv data/example_returns.csv 250 50 5 1 > $(BUILD_DIR)/fit_csv_output.txt
	$(BIN_DIR)/custom_constraint_example > $(BUILD_DIR)/constraint_output.txt

clean:
	rm -rf $(BUILD_DIR)
