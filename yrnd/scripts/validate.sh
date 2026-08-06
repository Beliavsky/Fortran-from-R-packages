#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
build="$root/build-validation"
rm -rf "$build"
mkdir -p "$build/debug" "$build/release"

compile_modules() {
  local out=$1
  shift
  local flags=("$@")
  cd "$out"
  local tvm="$root/dependencies/tvm-fortran/src"
  gfortran -std=f2018 -J. -I. "${flags[@]}" -c \
    "$tvm/tvm_kinds.f90" "$tvm/tvm_root.f90" \
    "$tvm/tvm_interpolation.f90" "$tvm/tvm_cashflows.f90" \
    "$tvm/tvm_curves.f90" "$tvm/tvm.f90" \
    "$root/src/yrnd_kinds.f90" "$root/src/yrnd_dates.f90" \
    "$root/src/yrnd_stats.f90" "$root/src/yrnd_optimize.f90" \
    "$root/src/yrnd_mixture.f90" "$root/src/yrnd_bonds.f90" \
    "$root/src/yrnd_transforms.f90" "$root/src/yrnd_api.f90" \
    "$root/src/yrnd.f90"
}

compile_and_run_tests() {
  local out=$1
  shift
  local flags=("$@")
  cd "$out"
  gfortran -std=f2018 -J. -I. "${flags[@]}" -c "$root/test/test_yrnd.f90"
  gfortran "${flags[@]}" -o test_yrnd *.o
  ./test_yrnd
  rm -f test_yrnd.o test_yrnd
  gfortran -std=f2018 -J. -I. "${flags[@]}" -c "$root/test/test_api.f90"
  gfortran "${flags[@]}" -o test_api *.o
  ./test_api
}

compile_modules "$build/debug" -O0 -g -Wall -Wextra -Wpedantic -fcheck=all -fbacktrace
compile_and_run_tests "$build/debug" -O0 -g -Wall -Wextra -Wpedantic -fcheck=all -fbacktrace

compile_modules "$build/release" -O3 -march=native -Wall -Wextra -Wpedantic
compile_and_run_tests "$build/release" -O3 -march=native -Wall -Wextra -Wpedantic

echo "yrnd validation completed successfully."
