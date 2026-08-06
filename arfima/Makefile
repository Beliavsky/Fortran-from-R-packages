FC = gfortran
SRC = src/arfima_kinds.f90 src/arfima_status.f90 src/arfima_types.f90 \
      src/arfima_linalg.f90 src/arfima_polynomial.f90 src/arfima_autocov.f90 \
      src/arfima_durbin.f90 src/arfima_transfer.f90 src/arfima_random.f90 \
      src/arfima_simulation.f90 src/arfima_forecast.f90 src/arfima_optimizer.f90 \
      src/arfima_information.f90 src/arfima_fit.f90 src/arfima_compat.f90 src/arfima.f90
TEST_NAMES = test_autocov test_durbin_sim test_forecast test_information_transfer test_fit test_extended
BASE = -std=f2018 -ffree-line-length-none -fimplicit-none -Wall -Wextra -Werror -Wno-compare-reals
CHECK = -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan
OPT = -O3

.PHONY: all check optimized example build-lib build-tests run-tests clean
all: check

check: FFLAGS=$(BASE) $(CHECK)
check: clean build-lib build-tests run-tests

optimized: FFLAGS=$(BASE) $(OPT)
optimized: clean build-lib build-tests run-tests

build-lib:
	@mkdir -p build/mod build/obj build/bin
	@cd build/obj && $(FC) $(FFLAGS) -J../mod -I../mod -c $(addprefix ../../,$(SRC))
	@$(FC) $(FFLAGS) -Jbuild/mod -Ibuild/mod -c test/test_support.f90 -o build/obj/test_support.o

build-tests:
	@for t in $(TEST_NAMES); do \
	  echo "Building $$t"; \
	  $(FC) $(FFLAGS) -Ibuild/mod build/obj/*.o test/$$t.f90 -o build/bin/$$t || exit 1; \
	done

run-tests:
	@for t in $(TEST_NAMES); do ./build/bin/$$t || exit 1; done

example: FFLAGS=$(BASE) $(OPT)
example: clean build-lib
	$(FC) $(FFLAGS) -Ibuild/mod build/obj/*.o example/arfima_demo.f90 -o build/bin/arfima_demo
	./build/bin/arfima_demo

clean:
	rm -rf build *.mod *.o src/*.mod src/*.o test/*.mod test/*.o example/*.mod example/*.o
