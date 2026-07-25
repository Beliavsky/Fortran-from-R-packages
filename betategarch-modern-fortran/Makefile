FC = gfortran
FFLAGS ?= -std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -ffree-line-length-none
CHECK_FLAGS ?= -std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror -ffree-line-length-none -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow

SRC = \
  src/betategarch_kinds.f90 \
  src/betategarch_math.f90 \
  src/betategarch_rng.f90 \
  src/skew_t_mod.f90 \
  src/bounded_nelder_mead_mod.f90 \
  src/tegarch_mod.f90 \
  src/betategarch.f90

BUILD = build
MODDIR = $(BUILD)/mod
OBJ = $(patsubst src/%.f90,$(BUILD)/%.o,$(SRC))

.PHONY: all check clean license-check

all: $(BUILD)/libbetategarch.a $(BUILD)/demo_betategarch $(BUILD)/fit_csv

$(BUILD) $(MODDIR):
	mkdir -p $@

$(BUILD)/%.o: src/%.f90 | $(BUILD) $(MODDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) -c $< -o $@

$(BUILD)/libbetategarch.a: $(OBJ)
	ar rcs $@ $^

$(BUILD)/demo_betategarch: app/demo_betategarch.f90 $(BUILD)/libbetategarch.a
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) $< $(BUILD)/libbetategarch.a -o $@

$(BUILD)/fit_csv: app/fit_csv.f90 $(BUILD)/libbetategarch.a
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) $< $(BUILD)/libbetategarch.a -o $@

$(BUILD)/example_simulate_and_fit: example/simulate_and_fit.f90 $(BUILD)/libbetategarch.a
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) $< $(BUILD)/libbetategarch.a -o $@

check: FFLAGS=$(CHECK_FLAGS)
check: clean license-check $(BUILD)/test_skew_t $(BUILD)/test_recursions $(BUILD)/test_models $(BUILD)/demo_betategarch $(BUILD)/fit_csv $(BUILD)/example_simulate_and_fit
	$(BUILD)/test_skew_t
	$(BUILD)/test_recursions
	$(BUILD)/test_models
	$(BUILD)/demo_betategarch > $(BUILD)/demo_output.txt
	$(BUILD)/fit_csv data/example_returns.csv 1 1 0 0 > $(BUILD)/fit_csv_output.txt
	$(BUILD)/example_simulate_and_fit > $(BUILD)/example_output.txt
	@echo "Application and example execution tests passed."

$(BUILD)/test_skew_t: test/test_skew_t.f90 $(BUILD)/libbetategarch.a
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) $< $(BUILD)/libbetategarch.a -o $@

$(BUILD)/test_recursions: test/test_recursions.f90 $(BUILD)/libbetategarch.a
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) $< $(BUILD)/libbetategarch.a -o $@

$(BUILD)/test_models: test/test_models.f90 $(BUILD)/libbetategarch.a
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) $< $(BUILD)/libbetategarch.a -o $@

license-check:
	sh test/check_license.sh

clean:
	rm -rf $(BUILD) *.mod
