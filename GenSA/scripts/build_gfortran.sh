#!/usr/bin/env bash
set -euo pipefail

mode=${1:-checked}
root=$(cd "$(dirname "$0")/.." && pwd)
build="$root/build-$mode"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"

common=(-std=f2018 -Wall -Wextra -Wpedantic -Wconversion-extra)
if [[ "$mode" == "checked" ]]; then
  flags=(-O0 -g -fcheck=all -fbacktrace)
elif [[ "$mode" == "optimized" ]]; then
  flags=(-O3)
else
  echo "usage: $0 [checked|optimized]" >&2
  exit 2
fi

sources=(
  gensa_kinds
  gensa_rng
  gensa_types
  gensa_local
  gensa
)
objects=()
for name in "${sources[@]}"; do
  object="$build/obj/$name.o"
  gfortran "${common[@]}" "${flags[@]}" -J"$build/mod" -I"$build/mod" \
    -c "$root/src/$name.f90" -o "$object"
  objects+=("$object")
done

for source in "$root"/test/*.f90; do
  name=$(basename "${source%.f90}")
  gfortran "${common[@]}" "${flags[@]}" -I"$build/mod" "$source" \
    "${objects[@]}" -o "$build/bin/$name"
  "$build/bin/$name"
done

for source in "$root"/example/*.f90 "$root"/app/*.f90; do
  name=$(basename "${source%.f90}")
  gfortran "${common[@]}" "${flags[@]}" -I"$build/mod" "$source" \
    "${objects[@]}" -o "$build/bin/$name"
  "$build/bin/$name" >/dev/null
done

echo "All $mode tests and examples passed."
