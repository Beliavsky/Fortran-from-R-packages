FC ?= gfortran
CHECK_FLAGS = -std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace -O0
OPT_FLAGS = -std=f2018 -Wall -Wextra -Werror -Wno-maybe-uninitialized -O3
BUILD_DIR = build

SRC = \
 src/rsdc_kinds.f90 src/rsdc_types.f90 src/rsdc_linalg.f90 \
 src/rsdc_random.f90 src/rsdc_parameters.f90 src/rsdc_filter.f90 \
 src/rsdc_likelihood.f90 \
 src/vendor/deoptimr_kinds.f90 src/vendor/deoptimr_interfaces.f90 \
 src/vendor/deoptimr_rng.f90 src/vendor/deoptimr_types.f90 \
 src/vendor/deoptimr_utils.f90 src/vendor/deoptimr_jde.f90 \
 src/vendor/deoptimr_ncde.f90 src/vendor/deoptimr.f90 \
 src/rsdc_simulation.f90 src/rsdc_forecast.f90 src/rsdc_viterbi.f90 \
 src/rsdc_estimation.f90 src/rsdc_portfolio.f90 src/rsdc_inference.f90 \
 src/rsdc_starts.f90 src/rsdc_bootstrap.f90 src/rsdc_bands.f90 src/rsdc.f90

TESTS = test_parameters test_filter test_simulation_forecast test_estimation \
 test_portfolio_inference test_statistics

.PHONY: check optimized demo clean

check:
	@rm -rf $(BUILD_DIR)/check && mkdir -p $(BUILD_DIR)/check
	@cd $(BUILD_DIR)/check && $(FC) $(CHECK_FLAGS) -c $(addprefix ../../,$(SRC))
	@for t in $(TESTS); do \
	  cd $(BUILD_DIR)/check && $(FC) $(CHECK_FLAGS) ../../test/$$t.f90 *.o -o $$t && ./$$t || exit 1; \
	  cd ../..; \
	done

optimized:
	@rm -rf $(BUILD_DIR)/optimized && mkdir -p $(BUILD_DIR)/optimized
	@cd $(BUILD_DIR)/optimized && $(FC) $(OPT_FLAGS) -c $(addprefix ../../,$(SRC))
	@for t in $(TESTS); do \
	  cd $(BUILD_DIR)/optimized && $(FC) $(OPT_FLAGS) ../../test/$$t.f90 *.o -o $$t && ./$$t || exit 1; \
	  cd ../..; \
	done

demo:
	@rm -rf $(BUILD_DIR)/demo && mkdir -p $(BUILD_DIR)/demo
	@cd $(BUILD_DIR)/demo && $(FC) $(OPT_FLAGS) -c $(addprefix ../../,$(SRC))
	@cd $(BUILD_DIR)/demo && $(FC) $(OPT_FLAGS) ../../app/rsdc_demo.f90 *.o -o rsdc_demo && ./rsdc_demo

clean:
	rm -rf $(BUILD_DIR)
