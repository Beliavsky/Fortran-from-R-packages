FC = gfortran
COMMON_FLAGS = -std=f2018 -Wall -Wextra -Wimplicit-interface -Werror -ffree-line-length-none
DEBUG_FLAGS = $(COMMON_FLAGS) -O0 -g -fcheck=all -fbacktrace
RELEASE_FLAGS = $(COMMON_FLAGS) -O2 -fbacktrace
LIBS = -llapack -lblas
SOURCES = \
 src/tsdyn_kinds.f90 \
 src/tsdyn_linalg.f90 \
 src/tsdyn_utils.f90 \
 src/tsdyn_ar.f90 \
 src/tsdyn_setar.f90 \
 src/tsdyn_lstar.f90 \
 src/tsdyn_local.f90 \
 src/tsdyn_var.f90 \
 src/tsdyn_vecm.f90 \
 src/tsdyn_tvar.f90 \
 src/tsdyn_tvecm.f90 \
 src/tsdyn_tests.f90 \
 src/tsdyn_girf.f90 \
 src/tsdyn_resample.f90 \
 src/tsdyn_rolling.f90 \
 src/tsdyn.f90
TESTS = test_univariate test_multivariate test_threshold
APPS = demo_tsdyn fit_csv threshold_example

.PHONY: all check debug-check release-check clean license
all: debug-check
license:
	./test/check_license.sh
check: debug-check

debug-check: license
	$(MAKE) clean-build
	$(MAKE) build-and-run MODE=debug FFLAGS='$(DEBUG_FLAGS)'

release-check: license
	$(MAKE) clean-build
	$(MAKE) build-and-run MODE=release FFLAGS='$(RELEASE_FLAGS)'

build-and-run:
	@mkdir -p build/$(MODE)/mod build/$(MODE)/obj build/$(MODE)/bin
	@set -e; for f in $(SOURCES); do \
	  b=$$(basename $$f .f90); \
	  $(FC) $(FFLAGS) -J build/$(MODE)/mod -I build/$(MODE)/mod -c $$f -o build/$(MODE)/obj/$$b.o; \
	done
	@$(FC) $(FFLAGS) -J build/$(MODE)/mod -I build/$(MODE)/mod -c test/test_support.f90 -o build/$(MODE)/obj/test_support.o
	@set -e; objs=$$(find build/$(MODE)/obj -name 'tsdyn_*.o' ! -name 'tsdyn.o' -print); \
	for t in $(TESTS); do \
	  $(FC) $(FFLAGS) -J build/$(MODE)/mod -I build/$(MODE)/mod -c test/$$t.f90 -o build/$(MODE)/obj/$$t.o; \
	  $(FC) $(FFLAGS) build/$(MODE)/obj/test_support.o build/$(MODE)/obj/$$t.o $$objs $(LIBS) -o build/$(MODE)/bin/$$t; \
	  build/$(MODE)/bin/$$t; \
	done
	@set -e; objs=$$(find build/$(MODE)/obj -name 'tsdyn_*.o' ! -name 'tsdyn.o' -print); \
	for a in demo_tsdyn fit_csv; do \
	  $(FC) $(FFLAGS) -J build/$(MODE)/mod -I build/$(MODE)/mod -c app/$$a.f90 -o build/$(MODE)/obj/$$a.o; \
	  $(FC) $(FFLAGS) build/$(MODE)/obj/$$a.o $$objs $(LIBS) -o build/$(MODE)/bin/$$a; \
	done; \
	$(FC) $(FFLAGS) -J build/$(MODE)/mod -I build/$(MODE)/mod -c example/threshold_example.f90 -o build/$(MODE)/obj/threshold_example.o; \
	$(FC) $(FFLAGS) build/$(MODE)/obj/threshold_example.o $$objs $(LIBS) -o build/$(MODE)/bin/threshold_example
	@build/$(MODE)/bin/demo_tsdyn >/dev/null
	@build/$(MODE)/bin/fit_csv data/sample.csv ar 2 >/dev/null
	@build/$(MODE)/bin/fit_csv data/sample.csv var 1 >/dev/null
	@build/$(MODE)/bin/threshold_example >/dev/null
	@echo '$(MODE) build, tests, and applications passed.'

clean-build:
	rm -rf build/debug build/release
clean:
	rm -rf build
