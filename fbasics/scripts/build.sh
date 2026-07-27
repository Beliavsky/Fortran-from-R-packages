#!/usr/bin/env sh
set -eu
mode=${1:-debug}
FC=${FC:-gfortran}
case "$mode" in
  debug) flags='-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror -Wno-unused-function -Wno-compare-reals -fcheck=all -fbacktrace -ffree-line-length-none' ;;
  release) flags='-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -Wno-unused-function -Wno-compare-reals -fbacktrace -ffree-line-length-none' ;;
  *) echo "unknown build mode: $mode" >&2; exit 2 ;;
esac
b="build/$mode"
rm -rf "$b"
mkdir -p "$b"
mods='src/fbasics_kinds.f90 src/fbasics_rng.f90 src/fbasics_special.f90 src/fbasics_linalg.f90 src/fbasics_stats.f90 src/fbasics_interp.f90 src/fbasics_optimize.f90 src/fbasics_distributions.f90 src/fbasics_stable.f90 src/fbasics_gld_extended.f90 src/fbasics_gh_fit.f90 src/fbasics_hypothesis.f90 src/fbasics_advanced_tests.f90 src/fbasics_moments.f90 src/fbasics_spatial.f90 src/fbasics_spline_density.f90 src/fbasics_drawdown.f90 src/fbasics_utils.f90 src/fbasics.f90 src/fbasics_test_support.f90'
for f in $mods; do
  o="$b/$(basename "${f%.f90}").o"
  "$FC" $flags -J"$b" -I"$b" -c "$f" -o "$o"
done
objs=''
for o in "$b"/fbasics_*.o; do objs="$objs $o"; done
for t in test/test_matrix_stats.f90 test/test_distributions.f90 test/test_inference_interp.f90 test/test_extended_algorithms.f90 test/test_gmm_prewhitening.f90; do
  exe="$b/$(basename "${t%.f90}")"
  "$FC" $flags -J"$b" -I"$b" "$t" "$b/fbasics.o" $objs -llapack -lblas -o "$exe"
done
for a in app/demo_fbasics.f90 app/analyze_csv.f90 example/distribution_example.f90; do
  exe="$b/$(basename "${a%.f90}")"
  "$FC" $flags -J"$b" -I"$b" "$a" "$b/fbasics.o" $objs -llapack -lblas -o "$exe"
done
