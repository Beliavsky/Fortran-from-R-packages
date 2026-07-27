FC = gfortran
FFLAGS ?= -O2 -std=f2018 -Wall -Wextra
CHECKFLAGS ?= -O0 -g -std=f2018 -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace
BUILD = build
SRC = src/rugarch_kinds.f90 src/rugarch_math.f90 src/rugarch_rng.f90 \
      src/rugarch_gh.f90 src/rugarch_distributions.f90 src/rugarch_optimizer.f90 \
      src/rugarch_linalg.f90 src/rugarch_resampling.f90 src/rugarch_types.f90 src/rugarch_models.f90 src/rugarch_fit.f90 \
      src/rugarch_risk.f90 src/rugarch_arfima.f90 \
      src/rugarch_backtests.f90 src/rugarch_inference.f90 src/rugarch_evaluation.f90 \
      src/rugarch_workflows.f90 src/rugarch_complete.f90 src/rugarch.f90

.PHONY: all test check demo clean

all: demo

demo:
	mkdir -p $(BUILD)
	cd $(BUILD) && $(FC) $(FFLAGS) $(addprefix ../,$(SRC)) ../app/demo_rugarch.f90 -o demo_rugarch

test:
	mkdir -p $(BUILD)
	cd $(BUILD) && for t in ../test/*.f90; do \
	  exe=$$(basename $$t .f90); \
	  $(FC) $(FFLAGS) $(addprefix ../,$(SRC)) $$t -o $$exe || exit 1; \
	  ./$$exe || exit 1; \
	done

check:
	mkdir -p $(BUILD)
	cd $(BUILD) && for t in ../test/*.f90; do \
	  exe=$$(basename $$t .f90); \
	  $(FC) $(CHECKFLAGS) $(addprefix ../,$(SRC)) $$t -o $$exe || exit 1; \
	  ./$$exe || exit 1; \
	done

clean:
	rm -rf $(BUILD) *.mod *.o
