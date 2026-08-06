#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
build="$root/build/example"
rm -rf "$build"
mkdir -p "$build"
cd "$build"

sources=(
  "$root/dependencies/splines/src/splines_kinds.f90"
  "$root/dependencies/splines/src/splines_linalg.f90"
  "$root/dependencies/splines/src/splines_core.f90"
  "$root/dependencies/splines/src/splines_basis.f90"
  "$root/dependencies/splines/src/splines.f90"
  "$root/src/mgcv_kinds.f90"
  "$root/src/mgcv_linalg.f90"
  "$root/src/mgcv_utils.f90"
  "$root/src/mgcv_distributions.f90"
  "$root/src/mgcv_smooths.f90"
  "$root/src/mgcv_families.f90"
  "$root/src/mgcv_fit.f90"
  "$root/src/mgcv_constraints.f90"
  "$root/src/mgcv_discrete.f90"
  "$root/src/mgcv_matrix.f90"
  "$root/src/mgcv_simulation.f90"
  "$root/src/mgcv.f90"
)

gfortran -std=f2018 -O2 -ffree-line-length-none \
  "${sources[@]}" "$root/example/basic_gam.f90" -o basic_gam
./basic_gam
