#!/usr/bin/env sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build-gfortran"
rm -rf "$build"
mkdir -p "$build"

flags="-std=f2018 -Wall -Wextra -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -ffree-line-length-none"
sources="
$root/src/riskportfolios_kinds.f90
$root/src/riskportfolios_linalg.f90
$root/src/riskportfolios_stats.f90
$root/src/riskportfolios_optimization.f90
$root/src/riskportfolios_mean.f90
$root/src/riskportfolios_semideviation.f90
$root/src/riskportfolios_covariance.f90
$root/src/riskportfolios_portfolio.f90
$root/src/riskportfolios.f90
"

gfortran $flags -J "$build" -I "$build" $sources \
  "$root/test/test_riskportfolios.f90" -llapack -lblas \
  -o "$build/test_riskportfolios"
"$build/test_riskportfolios"

gfortran $flags -J "$build" -I "$build" $sources \
  "$root/test/test_reference_values.f90" -llapack -lblas \
  -o "$build/test_reference_values"
"$build/test_reference_values"

gfortran -std=f2018 -O2 -Wall -Wextra -Werror -ffree-line-length-none \
  -J "$build" -I "$build" $sources "$root/app/demo_riskportfolios.f90" \
  -llapack -lblas -o "$build/demo_riskportfolios"
"$build/demo_riskportfolios"
