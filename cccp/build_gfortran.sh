#!/usr/bin/env bash
set -euo pipefail
mode=${1:-debug}
root=$(cd "$(dirname "$0")" && pwd)
build="$root/build-$mode"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
case "$mode" in
  debug) flags=(-O0 -g -std=f2018 -Wall -Wextra -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow) ;;
  release) flags=(-O3 -std=f2018 -Wall -Wextra -Werror) ;;
  *) echo "usage: $0 [debug|release]" >&2; exit 2 ;;
esac
sources=(cccp_kinds cccp_types cccp_linalg cccp_cones cccp_solver cccp_special cccp_problem cccp cccp_api)
objects=()
for name in "${sources[@]}"; do
  obj="$build/obj/$name.o"
  gfortran "${flags[@]}" -J "$build/mod" -I "$build/mod" -c "$root/src/$name.f90" -o "$obj"
  objects+=("$obj")
done
for src in "$root"/test/*.f90 "$root"/app/*.f90 "$root"/example/*.f90; do
  name=$(basename "${src%.f90}")
  gfortran "${flags[@]}" -J "$build/mod" -I "$build/mod" "$src" "${objects[@]}" -llapack -lblas -o "$build/bin/$name"
done
for exe in "$build"/bin/test_*; do "$exe"; done
"$build/bin/demo_cccp"
"$build/bin/second_order_cone"
"$build/bin/semidefinite_program"
