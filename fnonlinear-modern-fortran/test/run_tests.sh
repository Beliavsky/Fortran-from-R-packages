#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
mode=${1:-debug}
fc=${FC:-gfortran}
build="$root/build/$mode"
rm -rf "$build"
mkdir -p "$build"
common=(-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace)
case "$mode" in
  debug) flags=(-O0 -g -fcheck=all) ;;
  release) flags=(-O2) ;;
  *) echo "unknown mode: $mode" >&2; exit 2 ;;
esac
modflags=(-J"$build" -I"$build")
sources=(
  src/chaos_kinds.f90
  src/chaos_utils.f90
  src/chaos_embedding.f90
  src/chaos_metrics.f90
  src/chaos_neighbors.f90
  src/fnonlinear_rng.f90
  src/fnonlinear_maps.f90
  src/fnonlinear_linalg.f90
  src/fnonlinear_probability.f90
  src/fnonlinear_statistics.f90
  src/fnonlinear_tests.f90
  src/fnonlinear.f90
)
objects=()
for source in "${sources[@]}"; do
  object="$build/$(basename "${source%.f90}").o"
  "$fc" "${common[@]}" "${flags[@]}" "${modflags[@]}" -c "$root/$source" -o "$object"
  objects+=("$object")
done
for source in test/test_maps.f90 test/test_statistics.f90 test/test_neighbors.f90 test/test_tests.f90; do
  exe="$build/$(basename "${source%.f90}")"
  "$fc" "${common[@]}" "${flags[@]}" "${modflags[@]}" "$root/$source" "${objects[@]}" \
    -llapack -lblas -o "$exe"
  "$exe"
done
if [[ -f "$root/app/demo_fnonlinear.f90" ]]; then
  for source in app/demo_fnonlinear.f90 app/analyze_csv.f90 example/nonlinear_tests_example.f90; do
    exe="$build/$(basename "${source%.f90}")"
    "$fc" "${common[@]}" "${flags[@]}" "${modflags[@]}" "$root/$source" "${objects[@]}" \
      -llapack -lblas -o "$exe"
  done
  "$build/demo_fnonlinear" > "$build/demo.out"
  "$build/analyze_csv" "$root/data/logistic_sample.csv" > "$build/csv.out"
  "$build/nonlinear_tests_example" > "$build/example.out"
fi
"$root/test/check_license.sh"
echo "$mode build, tests, and applications passed."
