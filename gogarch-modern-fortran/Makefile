FC ?= gfortran

.PHONY: all check clean

all: check

check:
	FC="$(FC)" ./run_checks.sh

clean:
	rm -rf build
