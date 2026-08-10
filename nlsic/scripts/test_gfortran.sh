#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
build="$root/build-strict"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
flags=(-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all)
src=(nlsic_kinds nlsic_types nlsic_linalg nlsic_nnls nlsic_linear nlsic_solver nlsic)
for name in "${src[@]}"; do
  gfortran "${flags[@]}" -J "$build/mod" -I "$build/mod" \
    -c "$root/src/$name.f90" -o "$build/obj/$name.o"
done
objs=()
for name in "${src[@]}"; do objs+=("$build/obj/$name.o"); done
for file in "$root"/test/*.f90; do
  name=$(basename "${file%.f90}")
  gfortran "${flags[@]}" -I "$build/mod" "$file" "${objs[@]}" -o "$build/bin/$name"
  "$build/bin/$name"
done
for file in "$root"/example/*.f90; do
  name=$(basename "${file%.f90}")
  gfortran "${flags[@]}" -I "$build/mod" "$file" "${objs[@]}" -o "$build/bin/$name"
  "$build/bin/$name"
done
