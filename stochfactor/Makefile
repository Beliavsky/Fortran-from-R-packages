FC = gfortran
BASE = -std=f2018 -Wall -Wextra -Wimplicit-interface -Werror -ffree-line-length-none
DBG = -O0 -g -fcheck=all -fbacktrace
REL = -O2 -fbacktrace
LIBS = -llapack -lblas
SRC = src/sv_kinds.f90 src/sv_rng.f90 src/sv_linalg.f90 src/sv_stats.f90 src/sv_types.f90 src/sv_core.f90 src/fsv_types.f90 src/fsv_core.f90
BUILD = build
MOD = $(BUILD)/mod
.PHONY: all debug release test check clean apps
all: debug
debug: FLAGS=$(BASE) $(DBG)
debug: clean build_all test check apps
release: FLAGS=$(BASE) $(REL)
release: clean build_all test check apps
build_all:
	mkdir -p $(BUILD) $(MOD)
	$(FC) $(FLAGS) -J$(MOD) -I$(MOD) -c $(SRC)
	ar rcs $(BUILD)/libstochfactor.a *.o
	mv *.o $(BUILD)/
test: build_all
	$(FC) $(FLAGS) -J$(MOD) -I$(MOD) test/test_sv.f90 $(BUILD)/libstochfactor.a $(LIBS) -o $(BUILD)/test_sv
	$(FC) $(FLAGS) -J$(MOD) -I$(MOD) test/test_factor.f90 $(BUILD)/libstochfactor.a $(LIBS) -o $(BUILD)/test_factor
	$(BUILD)/test_sv
	$(BUILD)/test_factor
check:
	bash test/check_license.sh
apps: build_all
	$(FC) $(FLAGS) -J$(MOD) -I$(MOD) app/demo_sv.f90 $(BUILD)/libstochfactor.a $(LIBS) -o $(BUILD)/demo_sv
	$(FC) $(FLAGS) -J$(MOD) -I$(MOD) app/demo_factor.f90 $(BUILD)/libstochfactor.a $(LIBS) -o $(BUILD)/demo_factor
	$(FC) $(FLAGS) -J$(MOD) -I$(MOD) app/fit_csv.f90 $(BUILD)/libstochfactor.a $(LIBS) -o $(BUILD)/fit_csv
	$(FC) $(FLAGS) -J$(MOD) -I$(MOD) example/custom_restriction.f90 $(BUILD)/libstochfactor.a $(LIBS) -o $(BUILD)/custom_restriction
	$(BUILD)/demo_sv >/dev/null
	$(BUILD)/demo_factor >/dev/null
	$(BUILD)/fit_csv data/example.csv sv >/dev/null
	$(BUILD)/fit_csv data/example.csv factor 1 >/dev/null
	$(BUILD)/custom_restriction >/dev/null
clean:
	rm -rf $(BUILD) *.o *.mod
