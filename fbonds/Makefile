FC = gfortran
BUILD ?= build
MODE ?= debug
LIBS = -llapack -lblas
BASE_FLAGS = -std=f2018 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace
DEBUG_FLAGS = -O0 -g -fcheck=all
RELEASE_FLAGS = -O2
ifeq ($(MODE),release)
  FFLAGS = $(BASE_FLAGS) $(RELEASE_FLAGS)
else
  FFLAGS = $(BASE_FLAGS) $(DEBUG_FLAGS)
endif

MODULE_SOURCES = src/fbonds_kinds.f90 src/fbonds_linalg.f90 src/fbonds_csv.f90 src/fbonds_term_structure.f90
OBJECTS = $(BUILD)/fbonds_kinds.o $(BUILD)/fbonds_linalg.o $(BUILD)/fbonds_csv.o $(BUILD)/fbonds_term_structure.o

.PHONY: all test apps check debug release clean
all: apps test

debug:
	$(MAKE) clean
	$(MAKE) MODE=debug check

release:
	$(MAKE) clean
	$(MAKE) MODE=release check

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/fbonds_kinds.o: src/fbonds_kinds.f90 | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/fbonds_linalg.o: src/fbonds_linalg.f90 $(BUILD)/fbonds_kinds.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/fbonds_csv.o: src/fbonds_csv.f90 $(BUILD)/fbonds_kinds.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/fbonds_term_structure.o: src/fbonds_term_structure.f90 $(BUILD)/fbonds_kinds.o $(BUILD)/fbonds_linalg.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/test_term_structure: test/test_term_structure.f90 $(OBJECTS)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $^ $(LIBS) -o $@

$(BUILD)/demo_fbonds: app/demo_fbonds.f90 $(OBJECTS)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $^ $(LIBS) -o $@

$(BUILD)/fit_csv: app/fit_csv.f90 $(OBJECTS)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $^ $(LIBS) -o $@

$(BUILD)/fit_example: example/fit_example.f90 $(OBJECTS)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $^ $(LIBS) -o $@

test: $(BUILD)/test_term_structure
	./$(BUILD)/test_term_structure

apps: $(BUILD)/demo_fbonds $(BUILD)/fit_csv $(BUILD)/fit_example

check: apps test
	bash test/check_license.sh
	./$(BUILD)/demo_fbonds >/dev/null
	./$(BUILD)/fit_csv data/example_yield.csv ns >/dev/null
	./$(BUILD)/fit_csv data/example_yield.csv svensson l1 >/dev/null
	./$(BUILD)/fit_csv data/example_yield.csv svensson sse >/dev/null
	./$(BUILD)/fit_example >/dev/null
	@echo "$(MODE) build, tests, and applications passed."

clean:
	rm -rf $(BUILD) *.o *.mod
