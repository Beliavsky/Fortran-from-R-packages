FC ?= gfortran
BASE_FLAGS = -std=f2018 -ffree-line-length-none -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace
DEBUG_FLAGS = $(BASE_FLAGS) -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow
RELEASE_FLAGS = $(BASE_FLAGS) -O2

MODULES = kinds_mod probability_mod statistics_mod returns_mod drawdown_mod risk_mod capm_mod \
          performance_ratios_mod comoments_mod portfolio_risk_mod rolling_mod rng_mod \
          advanced_moments_mod tail_models_mod analytics_extensions_mod finite_sample_moments_mod csv_mod performanceanalytics_mod
TESTS = test_returns_drawdowns test_risk_ratios test_capm_moments test_rolling_cleaning test_advanced_algorithms test_finite_sample_corrections
APPS = demo_performanceanalytics analyze_csv portfolio_contributions advanced_estimators exact_moment_shrinkage

.PHONY: all debug release check clean
all: debug

debug:
	@./test/run_tests.sh debug

release:
	@./test/run_tests.sh release

check: debug
	@./test/check_license.sh

clean:
	rm -rf build
