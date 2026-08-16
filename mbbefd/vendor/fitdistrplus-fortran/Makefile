FC ?= gfortran
CHECK_FLAGS ?= -std=f2018 -O0 -Wall -Wextra -Werror -pedantic -fcheck=all -fbacktrace
OPT_FLAGS ?= -std=f2018 -O3 -Wall -Wextra -Werror -Wno-uninitialized -Wno-maybe-uninitialized -pedantic

SOURCES = \
  src/fitdistrplus_kinds.f90 \
  src/fitdistrplus_types.f90 \
  src/fitdistrplus_math.f90 \
  src/fitdistrplus_distributions.f90 \
  src/fitdistrplus_optimize.f90 \
  src/fitdistrplus_fit.f90 \
  src/fitdistrplus_stats.f90 \
  src/fitdistrplus_bootstrap.f90 \
  src/fitdistrplus.f90

TESTS = test_distributions test_fit_methods test_censored test_stats test_bootstrap_custom

.PHONY: all check optimized example clean
all: check

check:
	@$(MAKE) --no-print-directory run-suite BUILD=build/check FLAGS='$(CHECK_FLAGS)'

optimized:
	@$(MAKE) --no-print-directory run-suite BUILD=build/optimized FLAGS='$(OPT_FLAGS)'

run-suite:
	@rm -rf $(BUILD)
	@mkdir -p $(BUILD)
	@set -e; for source in $(SOURCES); do \
	  object=$(BUILD)/$$(basename $$source .f90).o; \
	  $(FC) $(FLAGS) -J $(BUILD) -I $(BUILD) -c $$source -o $$object; \
	done
	@$(FC) $(FLAGS) -J $(BUILD) -I $(BUILD) -c test/test_support.f90 -o $(BUILD)/test_support.o
	@set -e; for test_name in $(TESTS); do \
	  $(FC) $(FLAGS) -J $(BUILD) -I $(BUILD) $(BUILD)/*.o test/$$test_name.f90 -o $(BUILD)/$$test_name; \
	  $(BUILD)/$$test_name; \
	done

example:
	@rm -rf build/example
	@mkdir -p build/example
	@set -e; for source in $(SOURCES); do \
	  object=build/example/$$(basename $$source .f90).o; \
	  $(FC) $(OPT_FLAGS) -J build/example -I build/example -c $$source -o $$object; \
	done
	@$(FC) $(OPT_FLAGS) -J build/example -I build/example build/example/*.o \
	  app/demo_fitdistrplus.f90 -o build/example/demo_fitdistrplus
	@build/example/demo_fitdistrplus

clean:
	rm -rf build *.o *.mod
