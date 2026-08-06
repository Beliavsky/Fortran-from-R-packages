FC = gfortran
STD = -std=f2018
WARN = -Wall -Wextra -Werror -pedantic
CHECK_FLAGS = $(STD) $(WARN) -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
OPT_FLAGS = $(STD) $(WARN) -O3
SRC = src/vamc_kinds.f90 src/vamc_status.f90 src/vamc_dates.f90 src/vamc_math.f90 \
      src/vamc_curve.f90 src/vamc_scenarios.f90 src/vamc_policy.f90 \
      src/vamc_mortality.f90 src/vamc_valuation.f90 src/vamc.f90
TESTS = test_curve test_scenarios test_mortality_projection test_valuation test_portfolio_aging

.PHONY: check optimized example clean

check:
	@$(MAKE) MODE=check FLAGS='$(CHECK_FLAGS)' run-tests

optimized:
	@$(MAKE) MODE=optimized FLAGS='$(OPT_FLAGS)' run-tests

run-tests: build-objects
	$(FC) $(FLAGS) -J build/$(MODE) -I build/$(MODE) -c test/test_support.f90 -o build/$(MODE)/test_support.o
	@for t in $(TESTS); do \
	  $(FC) $(FLAGS) -J build/$(MODE) -I build/$(MODE) test/$$t.f90 build/$(MODE)/test_support.o build/$(MODE)/libvamc.a -o build/$(MODE)/$$t || exit 1; \
	  build/$(MODE)/$$t || exit 1; \
	done

build-objects:
	mkdir -p build/$(MODE)
	$(FC) $(FLAGS) -J build/$(MODE) -I build/$(MODE) -c $(SRC)
	mv *.o build/$(MODE)/
	ar rcs build/$(MODE)/libvamc.a build/$(MODE)/vamc_*.o build/$(MODE)/vamc.o

example: build-objects
	$(FC) $(FLAGS) -J build/$(MODE) -I build/$(MODE) app/vamc_demo.f90 build/$(MODE)/libvamc.a -o build/$(MODE)/vamc_demo
	build/$(MODE)/vamc_demo
	$(FC) $(FLAGS) -J build/$(MODE) -I build/$(MODE) example/curve_example.f90 build/$(MODE)/libvamc.a -o build/$(MODE)/curve_example
	build/$(MODE)/curve_example

clean:
	rm -rf build *.o *.mod
