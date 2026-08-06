FC = gfortran
SRC = src/kdensity_kinds.f90 src/kdensity_math.f90 src/kdensity_types.f90 \
      src/kdensity_starts.f90 src/kdensity_kernels.f90 src/kdensity_bandwidths.f90 \
      src/kdensity_core.f90 src/kdensity.f90
CHECK_FLAGS = -std=f2018 -O0 -g -Wall -Wextra -Werror -Wno-compare-reals -ftrampoline-impl=heap -ffree-line-length-none \
              -fcheck=all -fbacktrace -fimplicit-none -pedantic
OPT_FLAGS = -std=f2018 -O3 -Wall -Wextra -Werror -Wno-compare-reals -ftrampoline-impl=heap -ffree-line-length-none -fimplicit-none -pedantic
BUILD = build

.PHONY: check optimized example clean
check: clean
	mkdir -p $(BUILD)/check/mod $(BUILD)/check/obj $(BUILD)/check/bin
	$(FC) $(CHECK_FLAGS) -J$(BUILD)/check/mod -I$(BUILD)/check/mod -c $(SRC)
	mv *.o $(BUILD)/check/obj/
	for t in test/test_*.f90; do n=$$(basename $$t .f90); $(FC) $(CHECK_FLAGS) -I$(BUILD)/check/mod $$t $(BUILD)/check/obj/*.o -o $(BUILD)/check/bin/$$n; $(BUILD)/check/bin/$$n; done

optimized: clean
	mkdir -p $(BUILD)/opt/mod $(BUILD)/opt/obj $(BUILD)/opt/bin
	$(FC) $(OPT_FLAGS) -J$(BUILD)/opt/mod -I$(BUILD)/opt/mod -c $(SRC)
	mv *.o $(BUILD)/opt/obj/
	for t in test/test_*.f90; do n=$$(basename $$t .f90); $(FC) $(OPT_FLAGS) -I$(BUILD)/opt/mod $$t $(BUILD)/opt/obj/*.o -o $(BUILD)/opt/bin/$$n; $(BUILD)/opt/bin/$$n; done

example:
	mkdir -p $(BUILD)/example/mod $(BUILD)/example/obj
	$(FC) $(OPT_FLAGS) -J$(BUILD)/example/mod -I$(BUILD)/example/mod -c $(SRC)
	mv *.o $(BUILD)/example/obj/
	$(FC) $(OPT_FLAGS) -I$(BUILD)/example/mod example/kdensity_example.f90 $(BUILD)/example/obj/*.o -o $(BUILD)/example/kdensity_example
	$(BUILD)/example/kdensity_example

clean:
	rm -rf $(BUILD) *.o *.mod
