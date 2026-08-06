FC = gfortran
MODE ?= debug
BUILD_DIR := build/$(MODE)
MOD_DIR := $(BUILD_DIR)/mod
OBJ_DIR := $(BUILD_DIR)/obj
BIN_DIR := $(BUILD_DIR)/bin

SRC := src/matlab_kinds.f90 src/matlab_types.f90 src/matlab_numeric.f90 \
       src/matlab_array.f90 src/matlab_matrices.f90 src/matlab_number_theory.f90 \
       src/matlab_path.f90 src/matlab_timer.f90 src/matlab.f90
OBJ := $(patsubst src/%.f90,$(OBJ_DIR)/%.o,$(SRC))
TEST_SUPPORT := $(OBJ_DIR)/test_support.o
TESTS := test_numeric test_array test_matrices test_number_theory test_path_timer

ifeq ($(MODE),release)
  FFLAGS := -std=f2018 -O3 -Wall -Wextra -Werror -Wno-compare-reals -pedantic
else
  FFLAGS := -std=f2018 -O0 -g -Wall -Wextra -Werror -Wno-compare-reals -pedantic \
            -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow \
            -finit-real=snan -finit-integer=-999999
endif

.PHONY: all test example clean
all: test example

$(OBJ_DIR) $(MOD_DIR) $(BIN_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o: src/%.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(TEST_SUPPORT): test/test_support.f90 $(OBJ) | $(OBJ_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(BIN_DIR)/%: test/%.f90 $(OBJ) $(TEST_SUPPORT) | $(BIN_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(OBJ) $(TEST_SUPPORT) $< -o $@

$(BIN_DIR)/demo_matlab: example/demo_matlab.f90 $(OBJ) | $(BIN_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(OBJ) $< -o $@

test: $(addprefix $(BIN_DIR)/,$(TESTS))
	@set -e; for t in $(TESTS); do ./$(BIN_DIR)/$$t; done

example: $(BIN_DIR)/demo_matlab
	./$(BIN_DIR)/demo_matlab

clean:
	rm -rf build
