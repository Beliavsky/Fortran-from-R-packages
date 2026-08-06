FC = gfortran
BUILD ?= build

.PHONY: all check optimized demo clean

all: check

check:
	FC="$(FC)" BUILD="$(BUILD)/checked" ./scripts/test_gfortran.sh

optimized:
	FC="$(FC)" BUILD="$(BUILD)/optimized" ./scripts/test_gfortran_optimized.sh

demo:
	FC="$(FC)" BUILD="$(BUILD)/demo" ./scripts/run_demo.sh

clean:
	rm -rf $(BUILD) *.o *.mod
