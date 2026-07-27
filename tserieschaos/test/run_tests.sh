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
  src/chaos_systems.f90
  src/chaos_metrics.f90
  src/chaos_neighbors.f90
  src/tserieschaos.f90
)
objects=()
for source in "${sources[@]}"; do
  object="$build/$(basename "${source%.f90}").o"
  "$fc" "${common[@]}" "${flags[@]}" "${modflags[@]}" -c "$root/$source" -o "$object"
  objects+=("$object")
done
for source in test/test_core.f90 test/test_neighbors.f90 test/test_systems.f90; do
  exe="$build/$(basename "${source%.f90}")"
  "$fc" "${common[@]}" "${flags[@]}" "${modflags[@]}" "$root/$source" "${objects[@]}" -o "$exe"
  "$exe"
done
for source in app/demo_tserieschaos.f90 app/analyze_csv.f90 example/lorenz_analysis.f90; do
  exe="$build/$(basename "${source%.f90}")"
  "$fc" "${common[@]}" "${flags[@]}" "${modflags[@]}" "$root/$source" "${objects[@]}" -o "$exe"
done
"$build/demo_tserieschaos" > "$build/demo.out"
"$build/analyze_csv" "$root/data/lorenz_sample.csv" > "$build/csv_auto.out"
"$build/analyze_csv" "$root/data/lorenz_sample.csv" box > "$build/csv_box.out"
"$build/analyze_csv" "$root/data/lorenz_sample.csv" direct > "$build/csv_direct.out"
"$build/lorenz_analysis" > "$build/example.out"
"$root/test/check_license.sh"
echo "$mode build, tests, and applications passed."
