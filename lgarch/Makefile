FC = gfortran
BASE_FLAGS = -std=f2018 -Wall -Wextra -Wimplicit-interface -Werror -ffree-line-length-none
DEBUG_FLAGS = -O0 -g -fcheck=all -fbacktrace
RELEASE_FLAGS = -O2 -fbacktrace
LIBS = -llapack -lblas
SRC = src/lgarch_kinds.f90 src/lgarch_utils.f90 src/lgarch_rng.f90 src/lgarch_linalg.f90 src/lgarch_optimizer.f90 src/lgarch_univariate.f90 src/lgarch_multivariate.f90
BUILD = build
MODDIR = $(BUILD)/mod

.PHONY: all debug release check clean license apps
all: debug

debug: FFLAGS=$(BASE_FLAGS) $(DEBUG_FLAGS)
debug: $(BUILD)/liblgarch.a

release: clean
	$(MAKE) FFLAGS="$(BASE_FLAGS) $(RELEASE_FLAGS)" $(BUILD)/liblgarch.a

$(BUILD)/liblgarch.a: $(SRC)
	mkdir -p $(BUILD) $(MODDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) -c $(SRC)
	mv *.o $(BUILD)/
	ar rcs $@ $(BUILD)/*.o

license:
	bash test/check_license.sh

clean:
	rm -rf $(BUILD) *.mod *.o

TESTS = $(BUILD)/test_utils $(BUILD)/test_univariate $(BUILD)/test_multivariate
APPS = $(BUILD)/demo_lgarch $(BUILD)/fit_csv $(BUILD)/xreg_example

$(BUILD)/test_utils: test/test_utils.f90 $(BUILD)/liblgarch.a
	$(FC) $(FFLAGS) -I$(MODDIR) $< $(BUILD)/liblgarch.a $(LIBS) -o $@

$(BUILD)/test_univariate: test/test_univariate.f90 $(BUILD)/liblgarch.a
	$(FC) $(FFLAGS) -I$(MODDIR) $< $(BUILD)/liblgarch.a $(LIBS) -o $@

$(BUILD)/test_multivariate: test/test_multivariate.f90 $(BUILD)/liblgarch.a
	$(FC) $(FFLAGS) -I$(MODDIR) $< $(BUILD)/liblgarch.a $(LIBS) -o $@

check: FFLAGS=$(BASE_FLAGS) $(DEBUG_FLAGS)
check: license $(TESTS) $(APPS)
	$(BUILD)/test_utils
	$(BUILD)/test_univariate
	$(BUILD)/test_multivariate
	$(BUILD)/demo_lgarch > $(BUILD)/demo_output.txt
	$(BUILD)/xreg_example > $(BUILD)/xreg_output.txt
	$(BUILD)/fit_csv data/univariate.csv univariate ls 1 1 > $(BUILD)/csv_univariate_output.txt
	$(BUILD)/fit_csv data/multivariate.csv multivariate 0 0 > $(BUILD)/csv_multivariate_output.txt


$(BUILD)/demo_lgarch: app/demo_lgarch.f90 $(BUILD)/liblgarch.a
	$(FC) $(FFLAGS) -I$(MODDIR) $< $(BUILD)/liblgarch.a $(LIBS) -o $@

$(BUILD)/fit_csv: app/fit_csv.f90 $(BUILD)/liblgarch.a
	$(FC) $(FFLAGS) -I$(MODDIR) $< $(BUILD)/liblgarch.a $(LIBS) -o $@

$(BUILD)/xreg_example: example/xreg_example.f90 $(BUILD)/liblgarch.a
	$(FC) $(FFLAGS) -I$(MODDIR) $< $(BUILD)/liblgarch.a $(LIBS) -o $@

apps: $(APPS)
