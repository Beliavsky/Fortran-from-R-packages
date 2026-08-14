#!/usr/bin/env bash
set -euo pipefail
mode=${1:-checked}
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
case "$mode" in
  checked) flags=(-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -fbacktrace -O0) ;;
  optimized) flags=(-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -O3) ;;
  *) echo "usage: $0 checked|optimized" >&2; exit 2 ;;
esac
build="build/$mode"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
sources=(mco_kinds mco_random mco_pareto mco_quality mco_nsga2 mco_test_functions mco)
objects=()
for name in "${sources[@]}"; do
  gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" -c "src/$name.f90" -o "$build/obj/$name.o"
  objects+=("$build/obj/$name.o")
done
for source in test/*.f90 example/*.f90 app/*.f90; do
  name=$(basename "$source" .f90)
  gfortran "${flags[@]}" -I"$build/mod" "$source" "${objects[@]}" -o "$build/bin/$name"
  "$build/bin/$name"
done
