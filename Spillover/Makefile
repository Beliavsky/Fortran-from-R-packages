FC = gfortran
BUILD_DIR ?= build
MOD_DIR := $(BUILD_DIR)/mod
OBJ_DIR := $(BUILD_DIR)/obj
BIN_DIR := $(BUILD_DIR)/bin

BASE_FLAGS := -std=f2018 -Wall -Wextra -Werror
CHECK_FLAGS := $(BASE_FLAGS) -O0 -g -fcheck=all -fbacktrace
RELEASE_FLAGS := $(BASE_FLAGS) -O3
FFLAGS ?= $(CHECK_FLAGS)

SOURCES := src/spillover_kinds.f90 \
           src/spillover_status.f90 \
           src/spillover_linalg.f90 \
           src/spillover_random.f90 \
           src/spillover_var.f90 \
           src/spillover_fevd.f90 \
           src/spillover_tables.f90 \
           src/spillover_dynamic.f90 \
           src/spillover_compat.f90 \
           src/spillover.f90
OBJECTS := $(patsubst src/%.f90,$(OBJ_DIR)/%.o,$(SOURCES))
TEST_SOURCES := $(wildcard test/test_*.f90)
TEST_BINS := $(patsubst test/%.f90,$(BIN_DIR)/%,$(TEST_SOURCES))
EXAMPLE_BIN := $(BIN_DIR)/example_spillover

.PHONY: all check release test example clean directories

all: check

check: FFLAGS := $(CHECK_FLAGS)
check: clean test example

release: FFLAGS := $(RELEASE_FLAGS)
release: clean test example

test: $(TEST_BINS)
	@set -e; for test_program in $(TEST_BINS); do ./$$test_program; done

example: $(EXAMPLE_BIN)
	@./$(EXAMPLE_BIN)

$(OBJ_DIR)/spillover_kinds.o: src/spillover_kinds.f90 | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/spillover_status.o: src/spillover_status.f90 | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/spillover_linalg.o: src/spillover_linalg.f90 $(OBJ_DIR)/spillover_kinds.o $(OBJ_DIR)/spillover_status.o | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/spillover_random.o: src/spillover_random.f90 $(OBJ_DIR)/spillover_kinds.o | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/spillover_var.o: src/spillover_var.f90 $(OBJ_DIR)/spillover_kinds.o $(OBJ_DIR)/spillover_status.o $(OBJ_DIR)/spillover_linalg.o | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/spillover_fevd.o: src/spillover_fevd.f90 $(OBJ_DIR)/spillover_var.o | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/spillover_tables.o: src/spillover_tables.f90 $(OBJ_DIR)/spillover_fevd.o $(OBJ_DIR)/spillover_random.o | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/spillover_dynamic.o: src/spillover_dynamic.f90 $(OBJ_DIR)/spillover_tables.o | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/spillover_compat.o: src/spillover_compat.f90 $(OBJ_DIR)/spillover_dynamic.o | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/spillover.o: src/spillover.f90 $(OBJ_DIR)/spillover_compat.o | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(BIN_DIR)/test_%: test/test_%.f90 $(OBJECTS) | directories
	$(FC) $(FFLAGS) -I$(MOD_DIR) $< $(OBJECTS) -o $@

$(EXAMPLE_BIN): example/example_spillover.f90 $(OBJECTS) | directories
	$(FC) $(FFLAGS) -I$(MOD_DIR) $< $(OBJECTS) -o $@

directories:
	@mkdir -p $(MOD_DIR) $(OBJ_DIR) $(BIN_DIR)

clean:
	rm -rf $(BUILD_DIR)
