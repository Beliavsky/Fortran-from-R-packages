#!/usr/bin/env bash
set -euo pipefail
mode=${1:-strict}
root=$(cd "$(dirname "$0")" && pwd)
build="$root/build-gfortran-$mode"
rm -rf "$build" && mkdir -p "$build/mod" "$build/bin" "$build/obj"
common=(-std=f2018 -ffree-line-length-none -J "$build/mod" -I "$build/mod")
if [[ "$mode" == strict ]]; then
  flags=(-O0 -g -Wall -Wextra -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace)
elif [[ "$mode" == release ]]; then
  flags=(-O3 -Wall -Wextra -Werror)
else
  echo "usage: $0 [strict|release]" >&2; exit 2
fi
src=(
  src/bayesianou_kinds.f90 src/bayesianou_math.f90 src/bayesianou_types.f90
  src/bayesianou_utils.f90 src/bayesianou_geometry.f90 src/bayesianou_model.f90
  src/bayesianou_diagnostics.f90 src/bayesianou_mi.f90 src/bayesianou.f90
)
cd "$root"
objs=()
for file in "${src[@]}"; do
  obj="$build/obj/$(basename "${file%.f90}").o"
  gfortran "${common[@]}" "${flags[@]}" -c "$file" -o "$obj"
  objs+=("$obj")
done
for file in app/*.f90 example/*.f90 test/*.f90; do
  name=$(basename "${file%.f90}")
  gfortran "${common[@]}" "${flags[@]}" "$file" "${objs[@]}" -llapack -lblas -o "$build/bin/$name"
done
for exe in "$build"/bin/test_*; do "$exe"; done
"$build/bin/demo_bayesianou"
"$build/bin/nested_three_level"
"$build/bin/geometry_hmc"
echo "All $mode targets passed."
