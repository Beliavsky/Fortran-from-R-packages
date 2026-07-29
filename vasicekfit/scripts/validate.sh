#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
build="$root/build-validation"
rm -rf "$build"
mkdir -p "$build"
cd "$build"

sources=(
  "$root/src/vasicekfit_kinds.f90"
  "$root/src/vasicekfit_normal.f90"
  "$root/src/vasicekfit_linalg.f90"
  "$root/src/vasicekfit_distribution.f90"
  "$root/src/vasicekfit_model.f90"
  "$root/src/vasicekfit_inference.f90"
  "$root/src/vasicekfit.f90"
)
objects=(
  vasicekfit_kinds.o vasicekfit_normal.o vasicekfit_linalg.o
  vasicekfit_distribution.o vasicekfit_model.o vasicekfit_inference.o
  vasicekfit.o
)

run_build() {
  local optimization=$1
  shift
  local flags=(
    -std=f2018 -Wall -Wextra -Werror -Wconversion-extra
    -Wimplicit-interface "$optimization" "$@" -J . -I .
  )
  rm -f ./*.o ./*.mod ./*.exe ./test_* ./vasicekfit_demo \
    ./distribution_functions ./fit_and_predict
  gfortran "${flags[@]}" -c "${sources[@]}"
  for file in "$root"/test/*.f90 "$root"/app/*.f90 "$root"/example/*.f90; do
    name=$(basename "$file" .f90)
    gfortran "${flags[@]}" "$file" "${objects[@]}" -o "$name"
    "./$name"
  done
}

run_build -O0 -fcheck=all -fbacktrace
run_build -O2
printf '%s\n' 'validation: PASS'
