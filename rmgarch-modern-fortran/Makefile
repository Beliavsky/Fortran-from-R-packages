FC = gfortran
FFLAGS ?= -std=f2018 -O2 -Wall -Wextra -Wimplicit-interface
CHECKFLAGS ?= -std=f2018 -O0 -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace
BUILD := build
SRC := \
 src/rmgarch_kinds.f90 \
 src/rmgarch_types.f90 \
 src/rmgarch_math.f90 \
 src/rmgarch_rng.f90 \
 src/rmgarch_distributions.f90 \
 src/rmgarch_fft.f90 \
 src/rmgarch_optimizer.f90 \
 src/rmgarch_univariate.f90 \
 src/rmgarch_dcc.f90 \
 src/rmgarch_fdcc.f90 \
 src/rmgarch_scenario.f90 \
 src/rmgarch_copula.f90 \
 src/rmgarch_cgarch.f90 \
 src/rmgarch_ica.f90 \
 src/rmgarch_gogarch.f90 \
 src/rmgarch_varx.f90 \
 src/rmgarch_diagnostics.f90 \
 src/rmgarch_model.f90 \
 src/rmgarch_rolling.f90 \
 src/rmgarch.f90

.PHONY: all check check-license demo fit_csv clean

all: demo

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/librmgarch.a: $(SRC) | $(BUILD)
	cd $(BUILD) && $(FC) $(FFLAGS) -c $(addprefix ../,$(SRC))
	ar rcs $@ $(BUILD)/*.o

demo: $(BUILD)/librmgarch.a
	$(FC) $(FFLAGS) -I$(BUILD) app/demo_rmgarch.f90 $(BUILD)/librmgarch.a -o $(BUILD)/demo_rmgarch

fit_csv: $(BUILD)/librmgarch.a
	$(FC) $(FFLAGS) -I$(BUILD) example/fit_csv.f90 $(BUILD)/librmgarch.a -o $(BUILD)/fit_csv

check-license:
	sh test/check_license.sh

check: check-license | $(BUILD)
	cd $(BUILD) && $(FC) $(CHECKFLAGS) -c $(addprefix ../,$(SRC))
	$(FC) $(CHECKFLAGS) -I$(BUILD) test/test_dcc.f90 $(BUILD)/*.o -o $(BUILD)/test_dcc
	$(FC) $(CHECKFLAGS) -I$(BUILD) test/test_core.f90 $(BUILD)/*.o -o $(BUILD)/test_core
	$(FC) $(CHECKFLAGS) -I$(BUILD) test/test_extended.f90 $(BUILD)/*.o -o $(BUILD)/test_extended
	$(BUILD)/test_dcc
	$(BUILD)/test_core
	$(BUILD)/test_extended

clean:
	rm -rf $(BUILD) *.mod *.o
