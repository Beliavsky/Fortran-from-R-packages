FC = gfortran
FFLAGS ?= -O2 -std=f2018 -Wall -Wextra -Wpedantic
CHECK_FLAGS ?= -O0 -g -std=f2018 -Wall -Wextra -Wpedantic -fcheck=all -fbacktrace

SRC = src/longmemo_kinds.f90 \
      src/longmemo_fft.f90 \
      src/longmemo_linalg.f90 \
      src/longmemo_stats.f90 \
      src/longmemo_optimize.f90 \
      src/longmemo_io.f90 \
      src/longmemo.f90

BUILD = build
OBJECTS = $(BUILD)/longmemo_kinds.o \
          $(BUILD)/longmemo_fft.o \
          $(BUILD)/longmemo_linalg.o \
          $(BUILD)/longmemo_stats.o \
          $(BUILD)/longmemo_optimize.o \
          $(BUILD)/longmemo_io.o \
          $(BUILD)/longmemo.o

.PHONY: all test check strict-test example clean

all: $(BUILD)/liblongmemo.a

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/longmemo_kinds.o: src/longmemo_kinds.f90 | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/longmemo_fft.o: src/longmemo_fft.f90 $(BUILD)/longmemo_kinds.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/longmemo_linalg.o: src/longmemo_linalg.f90 $(BUILD)/longmemo_kinds.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/longmemo_stats.o: src/longmemo_stats.f90 $(BUILD)/longmemo_kinds.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/longmemo_optimize.o: src/longmemo_optimize.f90 $(BUILD)/longmemo_kinds.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/longmemo_io.o: src/longmemo_io.f90 $(BUILD)/longmemo_kinds.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/longmemo.o: src/longmemo.f90 $(BUILD)/longmemo_kinds.o $(BUILD)/longmemo_fft.o \
                      $(BUILD)/longmemo_linalg.o $(BUILD)/longmemo_stats.o $(BUILD)/longmemo_optimize.o $(BUILD)/longmemo_io.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/liblongmemo.a: $(OBJECTS)
	ar rcs $@ $(OBJECTS)

$(BUILD)/test_longmemo: $(OBJECTS) tests/test_longmemo.f90
	$(FC) $(CHECK_FLAGS) -J$(BUILD) -I$(BUILD) $(OBJECTS) tests/test_longmemo.f90 -o $@

test: $(BUILD)/test_longmemo
	./$(BUILD)/test_longmemo

check: strict-test

strict-test:
	rm -rf build-check
	mkdir -p build-check
	$(FC) $(CHECK_FLAGS) -Jbuild-check -Ibuild-check $(SRC) tests/test_longmemo.f90 -o build-check/test_longmemo
	./build-check/test_longmemo

$(BUILD)/demo_longmemo: $(OBJECTS) examples/demo_longmemo.f90
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $(OBJECTS) examples/demo_longmemo.f90 -o $@

example: $(BUILD)/demo_longmemo
	./$(BUILD)/demo_longmemo

clean:
	rm -rf $(BUILD) build-check *.mod *.o

$(BUILD)/fit_csv: $(OBJECTS) examples/fit_csv.f90
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $(OBJECTS) examples/fit_csv.f90 -o $@

fit-csv: $(BUILD)/fit_csv
	@echo "Run: ./$(BUILD)/fit_csv data/NileMin.csv"
