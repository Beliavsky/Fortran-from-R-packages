FC = gfortran
AR ?= ar
MODE ?= checked
BUILD := build/$(MODE)
MODDIR := $(BUILD)/mod
OBJDIR := $(BUILD)/obj
BINDIR := $(BUILD)/bin
LIB := $(BUILD)/libtvgarchkf.a

BASE_FLAGS := -std=f2018 -Wall -Wextra -Werror -pedantic
CHECKED_FLAGS := $(BASE_FLAGS) -O0 -g -fcheck=all -fbacktrace
OPTIMIZED_FLAGS := $(BASE_FLAGS) -O3
FFLAGS := $(if $(filter optimized,$(MODE)),$(OPTIMIZED_FLAGS),$(CHECKED_FLAGS))

SRC := src/fgarch_kinds.f90 src/fgarch_math.f90 src/fgarch_rng.f90 \
       src/fgarch_distributions.f90 src/fgarch_optimizer.f90 src/fgarch_types.f90 \
       src/fgarch_models.f90 src/fgarch_fit.f90 src/fgarch_risk.f90 src/fgarch.f90 \
       src/tvgarchkf_types.f90 src/tvgarchkf_functions.f90 \
       src/tvgarchkf_core.f90 src/tvgarchkf.f90
TESTS := test_functions test_filter test_simulation test_fit test_tv_parameter

.PHONY: all checked optimized test run-tests example clean

all: checked

$(LIB): $(SRC)
	mkdir -p $(MODDIR) $(OBJDIR) $(BINDIR)
	cd $(OBJDIR) && $(FC) $(FFLAGS) -J ../mod -I ../mod -c $(addprefix ../../../,$(SRC))
	$(AR) rcs $(LIB) $(OBJDIR)/*.o

$(OBJDIR)/test_support.o: test/manual/test_support.f90 $(LIB)
	$(FC) $(FFLAGS) -J $(MODDIR) -I $(MODDIR) -c $< -o $@

$(BINDIR)/%: test/manual/%.f90 $(OBJDIR)/test_support.o $(LIB)
	$(FC) $(FFLAGS) -J $(MODDIR) -I $(MODDIR) $< $(OBJDIR)/test_support.o $(LIB) -o $@

$(BINDIR)/demo_tvgarchkf: app/demo_tvgarchkf.f90 $(LIB)
	$(FC) $(FFLAGS) -J $(MODDIR) -I $(MODDIR) $< $(LIB) -o $@

run-tests: $(addprefix $(BINDIR)/,$(TESTS))
	@set -e; for test in $(TESTS); do echo "Running $$test"; $(BINDIR)/$$test; done

test: run-tests

example: $(BINDIR)/demo_tvgarchkf
	$(BINDIR)/demo_tvgarchkf

checked:
	$(MAKE) MODE=checked run-tests example

optimized:
	$(MAKE) MODE=optimized run-tests example

clean:
	rm -rf build *.o *.mod
