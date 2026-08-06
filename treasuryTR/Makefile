FC = gfortran
COMMON_FLAGS = -std=f2018 -Wall -Wextra -Wpedantic -Werror -Wno-compare-reals \
               -fimplicit-none -ffree-line-length-none
CHECK_FLAGS = $(COMMON_FLAGS) -O0 -g -fcheck=all -fbacktrace -finit-real=snan
OPT_FLAGS = $(COMMON_FLAGS) -O3

SOURCES = src/treasurytr_kinds.f90 \
          src/treasurytr_math.f90 \
          src/treasurytr_series.f90 \
          src/treasurytr.f90

TEST_NAMES = test_duration_convexity test_total_return test_matrix \
             test_preprocessing test_validation

.PHONY: all checked optimized test-checked test-optimized example clean

all: checked

checked:
	@mkdir -p build/checked/mod build/checked/bin
	@set -e; for name in $(TEST_NAMES); do \
	  $(FC) $(CHECK_FLAGS) -J build/checked/mod -I build/checked/mod \
	    $(SOURCES) test/$$name.f90 -o build/checked/bin/$$name; \
	done
	$(FC) $(CHECK_FLAGS) -J build/checked/mod -I build/checked/mod \
	  $(SOURCES) app/example_treasurytr.f90 -o build/checked/bin/example_treasurytr

optimized:
	@mkdir -p build/optimized/mod build/optimized/bin
	@set -e; for name in $(TEST_NAMES); do \
	  $(FC) $(OPT_FLAGS) -J build/optimized/mod -I build/optimized/mod \
	    $(SOURCES) test/$$name.f90 -o build/optimized/bin/$$name; \
	done
	$(FC) $(OPT_FLAGS) -J build/optimized/mod -I build/optimized/mod \
	  $(SOURCES) app/example_treasurytr.f90 -o build/optimized/bin/example_treasurytr

test-checked: checked
	@set -e; for name in $(TEST_NAMES); do ./build/checked/bin/$$name; done

test-optimized: optimized
	@set -e; for name in $(TEST_NAMES); do ./build/optimized/bin/$$name; done

example: checked
	./build/checked/bin/example_treasurytr

clean:
	rm -rf build *.mod *.o
