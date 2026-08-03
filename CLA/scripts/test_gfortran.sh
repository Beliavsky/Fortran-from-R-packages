#!/usr/bin/env bash
set -euo pipefail

mode=${1:-strict}
root=$(cd "$(dirname "$0")/.." && pwd)
build="$root/build-$mode"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"

case "$mode" in
  strict)
    flags=(-std=f2018 -O0 -g -Wall -Wextra -Werror -fcheck=all \
      -ffpe-trap=invalid,zero,overflow -fbacktrace)
    ;;
  optimized)
    flags=(-std=f2018 -O3 -Wall -Wextra -Werror)
    ;;
  *) echo "usage: $0 [strict|optimized]" >&2; exit 2 ;;
esac

sources=(
  kind_mod.f90
  cla_types.f90
  cla_core.f90
  cla_queries.f90
  cla_garch.f90
  cla.f90
)
objects=()
for src in "${sources[@]}"; do
  obj="$build/obj/${src%.f90}.o"
  gfortran "${flags[@]}" -J "$build/mod" -I "$build/mod" \
    -c "$root/src/$src" -o "$obj"
  objects+=("$obj")
done

run_target() {
  local file=$1
  local name
  name=$(basename "${file%.f90}")
  gfortran "${flags[@]}" -J "$build/mod" -I "$build/mod" \
    "$file" "${objects[@]}" -llapack -lblas -o "$build/bin/$name"
  "$build/bin/$name"
}

for file in "$root"/test/*.f90 "$root"/app/*.f90 "$root"/example/*.f90; do
  run_target "$file"
done

echo "CLA $mode build and tests passed"
