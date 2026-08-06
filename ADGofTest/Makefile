FC = gfortran
BUILD_DIR ?= build
MODE ?= release

BASE_FLAGS = -std=f2018 -pedantic -Wall -Wextra -Werror
CHECK_FLAGS = $(BASE_FLAGS) -O0 -g -fcheck=all -fbacktrace
RELEASE_FLAGS = $(BASE_FLAGS) -O3

ifeq ($(MODE),check)
  FFLAGS = $(CHECK_FLAGS)
else
  FFLAGS = $(RELEASE_FLAGS)
endif

MOD_DIR = $(BUILD_DIR)/$(MODE)/mod
OBJ_DIR = $(BUILD_DIR)/$(MODE)/obj
BIN_DIR = $(BUILD_DIR)/$(MODE)/bin
LIB_OBJ = $(OBJ_DIR)/adgoftest.o
TEST_SOURCES = $(wildcard test/*.f90)
TEST_BINS = $(patsubst test/%.f90,$(BIN_DIR)/%,$(TEST_SOURCES))
EXAMPLE_BIN = $(BIN_DIR)/adgoftest_demo

.PHONY: all test example check release clean

all: test example

check:
	$(MAKE) MODE=check all

release:
	$(MAKE) MODE=release all

test: $(TEST_BINS)
	@set -e; for test_program in $(TEST_BINS); do \
		echo "Running $$test_program"; \
		$$test_program; \
	done

example: $(EXAMPLE_BIN)
	$(EXAMPLE_BIN)

$(LIB_OBJ): src/adgoftest.f90
	@mkdir -p $(MOD_DIR) $(OBJ_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(BIN_DIR)/%: test/%.f90 $(LIB_OBJ)
	@mkdir -p $(BIN_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(LIB_OBJ) $< -o $@

$(EXAMPLE_BIN): example/adgoftest_demo.f90 $(LIB_OBJ)
	@mkdir -p $(BIN_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(LIB_OBJ) $< -o $@

clean:
	rm -rf $(BUILD_DIR)
