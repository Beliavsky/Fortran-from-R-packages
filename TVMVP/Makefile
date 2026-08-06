FC = gfortran
MODE ?= checked
BUILDDIR := build/$(MODE)
MODDIR := $(BUILDDIR)/mod
OBJDIR := $(BUILDDIR)/obj
BINDIR := $(BUILDDIR)/bin
SRC := \
 src/tvmvp_kinds.f90 \
 src/tvmvp_status.f90 \
 src/tvmvp_types.f90 \
 src/tvmvp_linalg.f90 \
 src/tvmvp_random.f90 \
 src/tvmvp_kernels.f90 \
 src/tvmvp_pca.f90 \
 src/tvmvp_poet.f90 \
 src/tvmvp_hypothesis.f90 \
 src/tvmvp_forecast.f90 \
 src/tvmvp_portfolio.f90 \
 src/tvmvp.f90
OBJ := $(patsubst src/%.f90,$(OBJDIR)/%.o,$(SRC))
BASEFLAGS := -std=f2018 -Wall -Wextra -Wpedantic -Werror -ffree-line-length-none -J$(MODDIR) -I$(MODDIR)
ifeq ($(MODE),checked)
  FFLAGS := $(BASEFLAGS) -O0 -g -fcheck=all -fbacktrace -finit-real=snan -finit-integer=-999999
else
  FFLAGS := $(BASEFLAGS) -O3
endif
TESTSRC := $(wildcard test/*.f90)
TESTBIN := $(patsubst test/%.f90,$(BINDIR)/%,$(TESTSRC))
EXAMPLESRC := $(wildcard example/*.f90)
EXAMPLEBIN := $(patsubst example/%.f90,$(BINDIR)/%,$(EXAMPLESRC))
.PHONY: all test example clean checked optimized
all: $(OBJ)
$(OBJDIR)/%.o: src/%.f90
	@mkdir -p $(OBJDIR) $(MODDIR) $(BINDIR)
	$(FC) $(FFLAGS) -c $< -o $@
$(BINDIR)/%: test/%.f90 $(OBJ)
	@mkdir -p $(BINDIR)
	$(FC) $(FFLAGS) $< $(OBJ) -o $@
$(BINDIR)/%: example/%.f90 $(OBJ)
	@mkdir -p $(BINDIR)
	$(FC) $(FFLAGS) $< $(OBJ) -o $@
test: $(TESTBIN)
	@set -e; for x in $(TESTBIN); do echo "Running $$x"; $$x; done
example: $(EXAMPLEBIN)
	@set -e; for x in $(EXAMPLEBIN); do echo "Running $$x"; $$x; done
checked:
	$(MAKE) MODE=checked test example
optimized:
	$(MAKE) MODE=optimized test example
clean:
	rm -rf build
