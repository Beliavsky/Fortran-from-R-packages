FC = gfortran
BUILD_DIR ?= build
MOD_DIR := $(BUILD_DIR)/mod
OBJ_DIR := $(BUILD_DIR)/obj
BIN_DIR := $(BUILD_DIR)/bin

BASE_FLAGS := -std=f2018 -Wall -Wextra -Werror
CHECK_FLAGS := $(BASE_FLAGS) -O0 -g -fcheck=all -fbacktrace
RELEASE_FLAGS := $(BASE_FLAGS) -O3
FFLAGS ?= $(CHECK_FLAGS)

SOURCES := src/smith_wilson_kinds.f90 \
           src/smith_wilson_linalg.f90 \
           src/smith_wilson.f90 \
           src/smith_wilson_yield_curve.f90
OBJECTS := $(OBJ_DIR)/smith_wilson_kinds.o \
           $(OBJ_DIR)/smith_wilson_linalg.o \
           $(OBJ_DIR)/smith_wilson.o \
           $(OBJ_DIR)/smith_wilson_yield_curve.o
TEST_SOURCES := $(wildcard test/test_*.f90)
TEST_BINS := $(patsubst test/%.f90,$(BIN_DIR)/%,$(TEST_SOURCES))
EXAMPLE_BIN := $(BIN_DIR)/example_smith_wilson

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

$(OBJ_DIR)/smith_wilson_kinds.o: src/smith_wilson_kinds.f90 | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/smith_wilson_linalg.o: src/smith_wilson_linalg.f90 $(OBJ_DIR)/smith_wilson_kinds.o | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/smith_wilson.o: src/smith_wilson.f90 $(OBJ_DIR)/smith_wilson_linalg.o | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/smith_wilson_yield_curve.o: src/smith_wilson_yield_curve.f90 $(OBJ_DIR)/smith_wilson.o | directories
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(BIN_DIR)/test_%: test/test_%.f90 $(OBJECTS) | directories
	$(FC) $(FFLAGS) -I$(MOD_DIR) $< $(OBJECTS) -o $@

$(EXAMPLE_BIN): example/example_smith_wilson.f90 $(OBJECTS) | directories
	$(FC) $(FFLAGS) -I$(MOD_DIR) $< $(OBJECTS) -o $@

directories:
	@mkdir -p $(MOD_DIR) $(OBJ_DIR) $(BIN_DIR)

clean:
	rm -rf $(BUILD_DIR)
