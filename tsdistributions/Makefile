FC = gfortran
BUILD ?= build
MODE ?= checked
BASE_FLAGS = -std=f2018 -pedantic -Wall -Wextra -Werror -Wconversion-extra -ffree-line-length-none
ifeq ($(MODE),checked)
  FFLAGS = $(BASE_FLAGS) -O0 -g -fcheck=all -fbacktrace
else
  FFLAGS = $(BASE_FLAGS) -O3
endif
MODFLAGS = -J$(BUILD) -I$(BUILD)

SOURCES = \
 src/ghyp_kinds.f90 src/ghyp_special.f90 src/ghyp_rng.f90 src/ghyp_linalg.f90 \
 src/ghyp_gig.f90 src/ghyp_model.f90 src/ghyp_distribution.f90 \
 src/tsd_types.f90 src/tsd_math.f90 src/tsd_distributions.f90 src/tsd_moments.f90 \
 src/tsd_optimize.f90 src/tsd_fit.f90 src/tsd_spd.f90 src/tsd_profile.f90 src/tsdistributions.f90
OBJECTS = $(patsubst src/%.f90,$(BUILD)/%.o,$(SOURCES))
TESTS = test_distributions test_random_moments test_estimation test_spd test_domain_profile
TEST_EXES = $(addprefix $(BUILD)/,$(TESTS))

.PHONY: all test example clean checked optimized
all: $(BUILD)/libtsdistributions.a

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/%.o: src/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAGS) -c $< -o $@

$(BUILD)/libtsdistributions.a: $(OBJECTS)
	ar rcs $@ $^

$(BUILD)/test_%: test/test_%.f90 $(BUILD)/libtsdistributions.a
	$(FC) $(FFLAGS) $(MODFLAGS) $< $(BUILD)/libtsdistributions.a -o $@

$(BUILD)/example_tsdistributions: example/example_tsdistributions.f90 $(BUILD)/libtsdistributions.a
	$(FC) $(FFLAGS) $(MODFLAGS) $< $(BUILD)/libtsdistributions.a -o $@

test: $(TEST_EXES)
	@set -e; for t in $(TEST_EXES); do ./$$t; done

example: $(BUILD)/example_tsdistributions
	./$(BUILD)/example_tsdistributions

checked:
	$(MAKE) clean
	$(MAKE) MODE=checked test example

optimized:
	$(MAKE) clean
	$(MAKE) MODE=optimized test example

clean:
	rm -rf $(BUILD)
