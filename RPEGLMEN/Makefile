FC = gfortran
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
  src/rpeglmen_kinds.f90 \
  src/rpeglmen_types.f90 \
  src/rpeglmen_math.f90 \
  src/rpeglmen_penalty.f90 \
  src/rpeglmen_likelihood.f90 \
  src/rpeglmen_solver.f90 \
  src/rpeglmen_gamma_mle.f90 \
  src/rpeglmen_cv.f90 \
  src/rpeglmen_compat.f90 \
  src/rpeglmen.f90
OBJECTS := $(patsubst src/%.f90,$(OBJ_DIR)/%.o,$(SOURCES))
TEST_NAMES := test_gradients test_exponential_fit test_gamma_fit test_cross_validation test_validation
TEST_BINS := $(addprefix $(BIN_DIR)/,$(TEST_NAMES))
EXAMPLE_BIN := $(BIN_DIR)/periodogram_glm_demo

.PHONY: all test example clean
all: $(TEST_BINS) $(EXAMPLE_BIN)

$(MOD_DIR) $(OBJ_DIR) $(BIN_DIR):
	mkdir -p $@

$(OBJ_DIR)/rpeglmen_kinds.o: src/rpeglmen_kinds.f90 | $(MOD_DIR) $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJ_DIR)/rpeglmen_types.o: src/rpeglmen_types.f90 $(OBJ_DIR)/rpeglmen_kinds.o
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJ_DIR)/rpeglmen_math.o: src/rpeglmen_math.f90 $(OBJ_DIR)/rpeglmen_kinds.o
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJ_DIR)/rpeglmen_penalty.o: src/rpeglmen_penalty.f90 $(OBJ_DIR)/rpeglmen_kinds.o
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJ_DIR)/rpeglmen_likelihood.o: src/rpeglmen_likelihood.f90 $(OBJ_DIR)/rpeglmen_kinds.o $(OBJ_DIR)/rpeglmen_math.o
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJ_DIR)/rpeglmen_solver.o: src/rpeglmen_solver.f90 $(OBJ_DIR)/rpeglmen_types.o $(OBJ_DIR)/rpeglmen_penalty.o $(OBJ_DIR)/rpeglmen_likelihood.o
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJ_DIR)/rpeglmen_gamma_mle.o: src/rpeglmen_gamma_mle.f90 $(OBJ_DIR)/rpeglmen_types.o $(OBJ_DIR)/rpeglmen_likelihood.o
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJ_DIR)/rpeglmen_cv.o: src/rpeglmen_cv.f90 $(OBJ_DIR)/rpeglmen_solver.o $(OBJ_DIR)/rpeglmen_gamma_mle.o
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJ_DIR)/rpeglmen_compat.o: src/rpeglmen_compat.f90 $(OBJ_DIR)/rpeglmen_cv.o
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJ_DIR)/rpeglmen.o: src/rpeglmen.f90 $(OBJ_DIR)/rpeglmen_compat.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BIN_DIR)/%: test/%.f90 $(OBJECTS) | $(BIN_DIR)
	$(FC) $(FFLAGS) $< $(OBJECTS) -o $@

$(EXAMPLE_BIN): example/periodogram_glm_demo.f90 $(OBJECTS) | $(BIN_DIR)
	$(FC) $(FFLAGS) $< $(OBJECTS) -o $@

test: $(TEST_BINS)
	@set -e; for test_program in $(TEST_BINS); do $$test_program; done

example: $(EXAMPLE_BIN)
	$(EXAMPLE_BIN)

clean:
	rm -rf build
