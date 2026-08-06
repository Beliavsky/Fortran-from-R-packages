FC = gfortran
BUILD ?= build/check
MODE ?= check

SOURCES = \
  src/vendor/rrcov_kinds.f90 \
  src/vendor/rrcov_types.f90 \
  src/vendor/rrcov_sort.f90 \
  src/vendor/rrcov_linalg.f90 \
  src/vendor/rrcov_stats.f90 \
  src/vendor/rrcov_random.f90 \
  src/vendor/rrcov_robust.f90 \
  src/vendor/rrcov_da.f90 \
  src/mass_types.f90 \
  src/mass_math.f90 \
  src/mass_basic.f90 \
  src/mass_regression.f90 \
  src/mass_robust.f90 \
  src/mass_discriminant.f90 \
  src/mass_distribution.f90 \
  src/mass_glm.f90 \
  src/mass_ordinal.f90 \
  src/mass_mds.f90 \
  src/mass_multivariate.f90 \
  src/mass_bandwidth.f90 \
  src/mass_model_selection.f90 \
  src/mass_nonlinear.f90 \
  src/mass.f90

ifeq ($(MODE),opt)
  FFLAGS = -std=f2018 -O3 -Wall -Wextra -Werror -pedantic \
    -Wno-uninitialized -Wno-maybe-uninitialized -J$(BUILD) -I$(BUILD)
else
  FFLAGS = -std=f2018 -O0 -g -Wall -Wextra -Werror -pedantic \
    -fcheck=all -fbacktrace -J$(BUILD) -I$(BUILD)
endif
LDFLAGS ?= -Wl,--no-warn-execstack

TESTS = test_basic test_regression test_distribution_glm \
  test_discriminant_robust test_multivariate_mds test_ordinal_bandwidth

.PHONY: all check optimized test example clean
all: check

check:
	$(MAKE) clean-build BUILD=build/check
	$(MAKE) test BUILD=build/check MODE=check

optimized:
	$(MAKE) clean-build BUILD=build/opt
	$(MAKE) test BUILD=build/opt MODE=opt

clean-build:
	rm -rf $(BUILD)
	mkdir -p $(BUILD)

$(BUILD)/modules.stamp: $(SOURCES)
	mkdir -p $(BUILD)
	@set -e; for source in $(SOURCES); do \
	  object=$(BUILD)/$$(basename $$source .f90).o; \
	  echo "Compiling $$source"; \
	  $(FC) $(FFLAGS) -c $$source -o $$object; \
	done
	touch $@

$(BUILD)/test_support.o: test/test_support.f90 $(BUILD)/modules.stamp
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/test_%: test/test_%.f90 $(BUILD)/test_support.o $(BUILD)/modules.stamp
	$(FC) $(FFLAGS) $< $(BUILD)/test_support.o \
	  $$(find $(BUILD) -maxdepth 1 -name '*.o' ! -name 'test_support.o' | sort) \
	  $(LDFLAGS) -o $@

test: $(addprefix $(BUILD)/,$(TESTS))
	@set -e; for test_program in $(addprefix $(BUILD)/,$(TESTS)); do \
	  echo "Running $$test_program"; $$test_program; \
	done

$(BUILD)/demo_mass: app/demo_mass.f90 $(BUILD)/modules.stamp
	$(FC) $(FFLAGS) $< $$(find $(BUILD) -maxdepth 1 -name '*.o' | sort) \
	  $(LDFLAGS) -o $@

example: $(BUILD)/demo_mass
	$(BUILD)/demo_mass

clean:
	rm -rf build *.mod *.o
