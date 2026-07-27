FC = gfortran
FFLAGS ?= -O2 -std=f2018 -Wall -Wextra -Wimplicit-interface
DEBUG_FLAGS ?= -O0 -g -std=f2018 -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace
BUILD = build
SRC = src/fgarch_kinds.f90 src/fgarch_math.f90 src/fgarch_rng.f90 \
      src/fgarch_distributions.f90 src/fgarch_optimizer.f90 src/fgarch_types.f90 \
      src/fgarch_models.f90 src/fgarch_fit.f90 src/fgarch_risk.f90 src/fgarch.f90

.PHONY: all test check clean

all: $(BUILD)/demo_fgarch

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/demo_fgarch: $(SRC) app/demo_fgarch.f90 | $(BUILD)
	$(FC) $(FFLAGS) -J $(BUILD) -I $(BUILD) $(SRC) app/demo_fgarch.f90 -o $@

$(BUILD)/test_distributions: $(SRC) test/test_distributions.f90 | $(BUILD)
	$(FC) $(DEBUG_FLAGS) -J $(BUILD) -I $(BUILD) $(SRC) test/test_distributions.f90 -o $@

$(BUILD)/test_models: $(SRC) test/test_models.f90 | $(BUILD)
	$(FC) $(DEBUG_FLAGS) -J $(BUILD) -I $(BUILD) $(SRC) test/test_models.f90 -o $@

test: $(BUILD)/test_distributions $(BUILD)/test_models
	$(BUILD)/test_distributions
	$(BUILD)/test_models

check: test $(BUILD)/demo_fgarch
	$(BUILD)/demo_fgarch

clean:
	rm -rf $(BUILD) *.mod
