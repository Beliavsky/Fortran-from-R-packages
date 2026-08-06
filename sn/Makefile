FC = gfortran
WARN_FLAGS = -std=f2018 -Wall -Wextra -Werror -fwrapv
CHECK_FLAGS = $(WARN_FLAGS) -O0 -g -fcheck=all -fbacktrace
OPT_FLAGS = $(WARN_FLAGS) -O3 -march=native
LINK_FLAGS = -Wl,--no-warn-execstack

MODULES = sn_kinds sn_status sn_math sn_rng sn_linalg sn_mvn \
          sn_univariate sn_multivariate sn_sun sn_matrix sn_parameters \
          sn_misc sn_optim sn_fit sn
SOURCES = $(addprefix src/,$(addsuffix .f90,$(MODULES)))
TESTS = test_fit test_multivariate test_rng_moments test_sun test_univariate test_utilities

.PHONY: all check optimized demo clean manifest

all: check

check:
	@rm -rf build/check && mkdir -p build/check
	@set -e; for f in $(SOURCES); do \
	  o=build/check/$$(basename $$f .f90).o; \
	  $(FC) $(CHECK_FLAGS) -J build/check -I build/check -c $$f -o $$o; \
	done
	@set -e; objs="$$(printf 'build/check/%s.o ' $(MODULES))"; \
	for t in $(TESTS); do \
	  $(FC) $(CHECK_FLAGS) -J build/check -I build/check test/$$t.f90 $$objs $(LINK_FLAGS) -o build/check/$$t; \
	  build/check/$$t; \
	done

optimized:
	@rm -rf build/optimized && mkdir -p build/optimized
	@set -e; for f in $(SOURCES); do \
	  o=build/optimized/$$(basename $$f .f90).o; \
	  $(FC) $(OPT_FLAGS) -J build/optimized -I build/optimized -c $$f -o $$o; \
	done
	@set -e; objs="$$(printf 'build/optimized/%s.o ' $(MODULES))"; \
	for t in $(TESTS); do \
	  $(FC) $(OPT_FLAGS) -J build/optimized -I build/optimized test/$$t.f90 $$objs $(LINK_FLAGS) -o build/optimized/$$t; \
	  build/optimized/$$t; \
	done

demo:
	@mkdir -p build/demo
	@set -e; for f in $(SOURCES); do \
	  o=build/demo/$$(basename $$f .f90).o; \
	  $(FC) $(OPT_FLAGS) -J build/demo -I build/demo -c $$f -o $$o; \
	done
	@objs="$$(printf 'build/demo/%s.o ' $(MODULES))"; \
	$(FC) $(OPT_FLAGS) -J build/demo -I build/demo app/demo_sn.f90 $$objs $(LINK_FLAGS) -o build/demo/demo_sn; \
	build/demo/demo_sn

manifest:
	@python3 -c "import tomllib; tomllib.load(open('fpm.toml','rb')); print('fpm.toml: valid TOML')"

clean:
	rm -rf build *.o *.mod
