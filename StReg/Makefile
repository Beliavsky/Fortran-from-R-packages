FC = gfortran
MODE ?= checked
BASE_FLAGS := -std=f2018 -ffree-line-length-none -Wall -Wextra -Wpedantic -Werror
CHECKED_FLAGS := -O0 -g -fcheck=all -fbacktrace
OPTIMIZED_FLAGS := -O3
ifeq ($(MODE),optimized)
  FFLAGS := $(BASE_FLAGS) $(OPTIMIZED_FLAGS)
else
  FFLAGS := $(BASE_FLAGS) $(CHECKED_FLAGS)
endif

SOURCES := src/streg_kinds.f90 src/streg_linalg.f90 src/streg_probability.f90 \
           src/streg_optimize.f90 src/streg_core.f90 src/streg.f90
TEST_DIRS := block_top conditional static dynamic inference_errors

.PHONY: all test example clean
all: test example

build:
	mkdir -p build

test: build
	@set -e; for d in $(TEST_DIRS); do \
	  echo "== test/$$d =="; \
	  $(FC) $(FFLAGS) -J build -I build $(SOURCES) test/$$d/test_support.f90 test/$$d/main.f90 -o build/test_$$d; \
	  ./build/test_$$d; \
	done

example: build
	$(FC) $(FFLAGS) -J build -I build $(SOURCES) example/streg_demo.f90 -o build/streg_demo
	./build/streg_demo

clean:
	rm -rf build *.mod *.o
