#!/usr/bin/env sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${TMPDIR:-/tmp}/mts-fortran-validate-$$
trap 'rm -rf "$build"' EXIT INT TERM
mkdir -p "$build"
cd "$build"

fc=${FC:-gfortran}
flags='-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace -ffree-line-length-none'
modules='mts_kinds mts_types mts_linalg mts_stats mts_rng mts_optimize mts_var mts_varma mts_regression mts_diagnostics mts_volatility mts_factor mts_ecm_missing mts_structural mts'
objects=''
for module in $modules; do
    "$fc" $flags -J . -I . -c "$root/src/$module.f90"
    objects="$objects $module.o"
done
"$fc" $flags -J . -I . -c "$root/test/test_support.f90"

for test in test_linalg_structural test_var test_varma_varx test_diagnostics_factor test_volatility test_ecm_missing; do
    "$fc" $flags -J . -I . "$root/test/$test.f90" test_support.o $objects -o "$test"
    "./$test"
done

"$fc" $flags -J . -I . "$root/app/demo_mts.f90" $objects -o demo_mts
./demo_mts >/dev/null
for example in var_analysis varma_analysis volatility_analysis factor_vecm_analysis structural_tools; do
    "$fc" $flags -J . -I . "$root/example/$example.f90" $objects -o "$example"
    "./$example" >/dev/null
done

echo 'checked validation: PASS'
