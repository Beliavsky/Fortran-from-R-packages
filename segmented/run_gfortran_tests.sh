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
  nlme_lme nlme_nonlinear nlme_diagnostics nlme_grouped nlme_models nlme
  segmented_status segmented_types segmented_utils segmented_fit
  segmented_mixed segmented_inference segmented_wrappers segmented
)
pushd "$build/obj" >/dev/null
for source in "${sources[@]}"; do
  gfortran "${flags[@]}" -c "$root/src/$source.f90" -J "$build/mod" -I "$build/mod"
done
for test in test_segmented_lm test_stepmented_glm test_inference_selection test_segmented_lme; do
  gfortran "${flags[@]}" "$root/test/$test.f90" ./*.o \
    -J "$build/mod" -I "$build/mod" -o "$build/bin/$test"
  "$build/bin/$test"
done
for example in segmented_linear_example stepmented_example segmented_glm_example segmented_lme_example; do
  gfortran "${flags[@]}" "$root/example/$example.f90" ./*.o \
    -J "$build/mod" -I "$build/mod" -o "$build/bin/$example"
  "$build/bin/$example" >/dev/null
done
gfortran "${flags[@]}" "$root/app/demo_segmented.f90" ./*.o \
  -J "$build/mod" -I "$build/mod" -o "$build/bin/demo_segmented"
"$build/bin/demo_segmented" >/dev/null
popd >/dev/null
echo "GNU Fortran $mode validation: PASS"
