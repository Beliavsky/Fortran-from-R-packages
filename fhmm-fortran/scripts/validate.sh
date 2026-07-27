#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")/.."
rm -rf build_validation
mkdir build_validation
flags="-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace"
sources="src/fhmm_kinds.f90 src/fhmm_types.f90 src/fhmm_math.f90 src/fhmm_distributions.f90 src/fhmm_parameters.f90 src/fhmm_algorithms.f90 src/fhmm_hierarchical.f90 src/fhmm_optimize.f90 src/fhmm_estimation.f90 src/fhmm_diagnostics.f90 src/fhmm_calendar.f90 src/fhmm.f90"
gfortran $flags -J build_validation -I build_validation -c $sources
for file in test/*.f90; do
  exe="build_validation/$(basename "${file%.f90}")"
  gfortran $flags -J build_validation -I build_validation ./*.o "$file" -o "$exe"
  "$exe"
done
for file in app/*.f90 example/*.f90; do
  exe="build_validation/$(basename "${file%.f90}")"
  gfortran $flags -J build_validation -I build_validation ./*.o "$file" -o "$exe"
  "$exe" >/dev/null
done
rm -f ./*.o ./*.mod
rm -rf build_validation
printf '%s\n' 'validation: PASS'
