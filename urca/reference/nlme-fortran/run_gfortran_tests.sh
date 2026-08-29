#!/usr/bin/env bash
set -euo pipefail
mode="${1:-debug}"
root="$(cd "$(dirname "$0")" && pwd)"
build="$root/build-gfortran-$mode"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
common=(-std=f2018 -Wall -Wextra -Wpedantic -Werror -fbacktrace)
if [[ "$mode" == "release" ]]; then
  flags=("${common[@]}" -O3)
else
  flags=("${common[@]}" -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow)
fi
sources=(
  nlme_kinds nlme_status nlme_types nlme_linalg nlme_correlation
  nlme_variance nlme_pdmat nlme_optimize nlme_covariance nlme_gls
  nlme_lme nlme_nonlinear nlme_diagnostics nlme_grouped nlme_models
  nlme_test_support nlme
)
pushd "$build/obj" >/dev/null
for source in "${sources[@]}"; do
  gfortran "${flags[@]}" -c "$root/src/$source.f90" -J "$build/mod" -I "$build/mod"
done
for test in test_covariance test_gls test_lme test_nonlinear_diagnostics; do
  gfortran "${flags[@]}" "$root/test/$test.f90" ./*.o \
    -J "$build/mod" -I "$build/mod" -o "$build/bin/$test"
  "$build/bin/$test"
done
for example in gls_ar1_example random_intercept_lme nonlinear_gnls_example correlation_matrices; do
  gfortran "${flags[@]}" "$root/example/$example.f90" ./*.o \
    -J "$build/mod" -I "$build/mod" -o "$build/bin/$example"
  "$build/bin/$example" >/dev/null
done
gfortran "${flags[@]}" "$root/app/demo_nlme.f90" ./*.o \
  -J "$build/mod" -I "$build/mod" -o "$build/bin/demo_nlme"
"$build/bin/demo_nlme" >/dev/null
popd >/dev/null
echo "GNU Fortran $mode validation: PASS"
