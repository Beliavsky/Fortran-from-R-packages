FC = gfortran
SRC = src/ltsa_kinds.f90 src/ltsa_status.f90 src/ltsa_types.f90 src/ltsa_random.f90 \
      src/ltsa_linalg.f90 src/ltsa_fft.f90 src/ltsa_durbin_levinson.f90 \
      src/ltsa_toeplitz.f90 src/ltsa_arma.f90 src/ltsa_simulation.f90 \
      src/ltsa_likelihood.f90 src/ltsa_forecast.f90 src/ltsa_innovation.f90 src/ltsa_compat.f90 src/ltsa.f90
TESTS = test/test_durbin test/test_toeplitz test/test_arma_forecast test/test_simulation test/test_innovation
BASE = -std=f2018 -ffree-line-length-none -fimplicit-none -Wall -Wextra -Werror -Wno-compare-reals
CHECK = -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan
OPT = -O3 -march=native

.PHONY: all check optimized example clean
all: check

check: FFLAGS=$(BASE) $(CHECK)
check: clean build-tests run-tests

optimized: FFLAGS=$(BASE) $(OPT)
optimized: clean build-tests run-tests

build-tests:
	@mkdir -p build/mod build/bin
	$(FC) $(FFLAGS) -Jbuild/mod -Ibuild/mod $(SRC) test/test_durbin.f90 -o build/bin/test_durbin
	$(FC) $(FFLAGS) -Jbuild/mod -Ibuild/mod $(SRC) test/test_toeplitz.f90 -o build/bin/test_toeplitz
	$(FC) $(FFLAGS) -Jbuild/mod -Ibuild/mod $(SRC) test/test_arma_forecast.f90 -o build/bin/test_arma_forecast
	$(FC) $(FFLAGS) -Jbuild/mod -Ibuild/mod $(SRC) test/test_simulation.f90 -o build/bin/test_simulation
	$(FC) $(FFLAGS) -Jbuild/mod -Ibuild/mod $(SRC) test/test_innovation.f90 -o build/bin/test_innovation

run-tests:
	@for t in $(TESTS); do ./build/bin/$$(basename $$t); done

example: FFLAGS=$(BASE) $(OPT)
example: clean
	@mkdir -p build/mod build/bin
	$(FC) $(FFLAGS) -Jbuild/mod -Ibuild/mod $(SRC) example/ltsa_demo.f90 -o build/bin/ltsa_demo
	./build/bin/ltsa_demo

clean:
	rm -rf build *.mod src/*.mod test/*.mod example/*.mod
