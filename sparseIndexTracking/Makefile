FC = gfortran
BUILD_DIR ?= build
MOD_DIR := $(BUILD_DIR)/mod
OBJ_DIR := $(BUILD_DIR)/obj
BIN_DIR := $(BUILD_DIR)/bin

BASE_FLAGS := -std=f2018 -Wall -Wextra -Werror
CHECK_FLAGS := $(BASE_FLAGS) -O0 -g -fcheck=all -fbacktrace
RELEASE_FLAGS := $(BASE_FLAGS) -O3
FFLAGS ?= $(CHECK_FLAGS)

SOURCES := src/sparse_index_tracking_kinds.f90 \
           src/sparse_index_tracking_linalg.f90 \
           src/sparse_index_tracking_projection.f90 \
           src/sparse_index_tracking_core.f90 \
           src/sparse_index_tracking.f90
OBJECTS := $(OBJ_DIR)/sparse_index_tracking_kinds.o \
           $(OBJ_DIR)/sparse_index_tracking_linalg.o \
           $(OBJ_DIR)/sparse_index_tracking_projection.o \
           $(OBJ_DIR)/sparse_index_tracking_core.o \
           $(OBJ_DIR)/sparse_index_tracking.o
TEST_SOURCES := $(wildcard test/test_*.f90)
TEST_BINS := $(patsubst test/%.f90,$(BIN_DIR)/%,$(TEST_SOURCES))
EXAMPLE_BIN := $(BIN_DIR)/example_sparse_index_tracking

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

$(OBJ_DIR)/sparse_index_tracking_kinds.o: src/sparse_index_tracking_kinds.f90 | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/sparse_index_tracking_linalg.o: src/sparse_index_tracking_linalg.f90 $(OBJ_DIR)/sparse_index_tracking_kinds.o | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/sparse_index_tracking_projection.o: src/sparse_index_tracking_projection.f90 $(OBJ_DIR)/sparse_index_tracking_linalg.o | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/sparse_index_tracking_core.o: src/sparse_index_tracking_core.f90 $(OBJ_DIR)/sparse_index_tracking_projection.o | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/sparse_index_tracking.o: src/sparse_index_tracking.f90 $(OBJ_DIR)/sparse_index_tracking_core.o | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(BIN_DIR)/test_%: test/test_%.f90 $(OBJECTS) | directories
	$(FC) $(FFLAGS) -I$(MOD_DIR) $< $(OBJECTS) -o $@

$(EXAMPLE_BIN): example/example_sparse_index_tracking.f90 $(OBJECTS) | directories
	$(FC) $(FFLAGS) -I$(MOD_DIR) $< $(OBJECTS) -o $@

directories:
	@mkdir -p $(MOD_DIR) $(OBJ_DIR) $(BIN_DIR)

clean:
	rm -rf $(BUILD_DIR)
