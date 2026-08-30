#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
work=${TMPDIR:-/tmp}/gbm3-fortran-strict-$$
mkdir -p "$work/mod"
trap 'rm -rf "$work"' EXIT INT TERM
cd "$root"
flags='-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fimplicit-none -fcheck=all'
gfortran -c $flags -J "$work/mod" -I "$work/mod" \
  src/gbm3_kinds.f90 \
  src/gbm3_constants.f90 \
  src/gbm3_math.f90 \
  src/gbm3_types.f90 \
  src/gbm3_tree.f90 \
  src/gbm3_distributions.f90 \
  src/gbm3_pairwise.f90 \
  src/gbm3_cox.f90 \
  src/gbm3_core.f90 \
  src/gbm3_cv.f90 \
  src/gbm3_diagnostics.f90 \
  src/gbm3.f90
gfortran $flags -I "$work/mod" -J "$work/mod" test/test_gbm3.f90 ./*.o -o "$work/test_gbm3"
"$work/test_gbm3"
gfortran $flags -I "$work/mod" -J "$work/mod" example/gaussian_example.f90 ./*.o -o "$work/gaussian_example"
"$work/gaussian_example"
gfortran $flags -I "$work/mod" -J "$work/mod" example/advanced_api_example.f90 ./*.o -o "$work/advanced_api_example"
"$work/advanced_api_example"
python3 tools/check_source.py
rm -f ./*.o
