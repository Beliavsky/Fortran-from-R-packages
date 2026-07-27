FC = gfortran
FFLAGS ?= -std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace
SRC = src/tseries_kinds.f90 src/tseries_types.f90 src/tseries_special.f90 \
      src/tseries_stats.f90 src/tseries_linalg.f90 src/tseries_random.f90 \
      src/tseries_fft.f90 src/tseries_optimize.f90 src/tseries_models.f90 \
      src/tseries_tests.f90 src/tseries_resampling.f90 src/tseries_finance.f90 src/tseries.f90

.PHONY: all test clean
all: demo_tseries

demo_tseries: $(SRC) app/demo_tseries.f90 | build
	$(FC) $(FFLAGS) -J build -I build $(SRC) app/demo_tseries.f90 -o $@

test_smoke: $(SRC) test/test_smoke.f90 | build
	$(FC) $(FFLAGS) -J build -I build $(SRC) test/test_smoke.f90 -o $@

test: test_smoke test_extended
	./test_smoke
	./test_extended

test_extended: $(SRC) test/test_extended.f90 | build
	$(FC) $(FFLAGS) -J build -I build $(SRC) test/test_extended.f90 -o $@

build:
	mkdir -p build

clean:
	rm -rf build *.o *.mod demo_tseries test_smoke test_extended
