FC = gfortran
AR = ar
BUILD ?= build
MODDIR := $(BUILD)/mod
OBJDIR := $(BUILD)/obj
LIB := $(BUILD)/libpolynom.a

SOURCES := \
 src/polynom_kinds.f90 \
 src/polynom_status.f90 \
 src/polynom_core.f90 \
 src/polynom_algorithms.f90 \
 src/polynom_roots.f90 \
 src/polynom_orthogonal.f90 \
 src/polynom_compat.f90 \
 src/polynom.f90
OBJECTS := $(patsubst src/%.f90,$(OBJDIR)/%.o,$(SOURCES))
TESTS := test_arithmetic test_calculus test_orthogonal test_gcd_summary test_api_edges

BASE_FLAGS := -std=f2018 -Wall -Wextra -Werror -pedantic
CHECK_FLAGS := $(BASE_FLAGS) -O0 -g -fcheck=all -fbacktrace
OPT_FLAGS := $(BASE_FLAGS) -O3
FFLAGS ?= $(CHECK_FLAGS)

.PHONY: all checked optimized test example clean
all: checked

$(OBJDIR)/%.o: src/%.f90
	@mkdir -p $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) -c $< -o $@

$(LIB): $(OBJECTS)
	@mkdir -p $(BUILD)
	$(AR) rcs $@ $^

$(BUILD)/%: test/%.f90 $(LIB)
	$(FC) $(FFLAGS) -I$(MODDIR) $< $(LIB) -o $@

$(BUILD)/demo_polynom: example/demo_polynom.f90 $(LIB)
	$(FC) $(FFLAGS) -I$(MODDIR) $< $(LIB) -o $@

checked:
	$(MAKE) clean
	$(MAKE) FFLAGS='$(CHECK_FLAGS)' test example

optimized:
	$(MAKE) clean
	$(MAKE) FFLAGS='$(OPT_FLAGS)' test example

test: $(addprefix $(BUILD)/,$(TESTS))
	@set -e; for t in $(TESTS); do ./$(BUILD)/$$t; done

example: $(BUILD)/demo_polynom
	./$(BUILD)/demo_polynom

clean:
	rm -rf build build_check build_extract
