#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
rm -rf build-strict
mkdir build-strict
flags=(-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0)
mods=(-Jbuild-strict -Ibuild-strict)
srcs=(src/lpsolve_types.f90 src/lpsolve_simplex.f90 src/lpsolve_core.f90 src/lpsolve_special.f90 src/lpsolve.f90)
objs=()
for src in "${srcs[@]}"; do
  obj="build-strict/$(basename "${src%.f90}").o"
  gfortran "${flags[@]}" "${mods[@]}" -c "$src" -o "$obj"
  objs+=("$obj")
done
for test in test/*.f90; do
  exe="build-strict/$(basename "${test%.f90}")"
  gfortran "${flags[@]}" "${mods[@]}" "$test" "${objs[@]}" -o "$exe"
  "$exe"
done
for ex in example/*.f90; do
  exe="build-strict/$(basename "${ex%.f90}")"
  gfortran "${flags[@]}" "${mods[@]}" "$ex" "${objs[@]}" -o "$exe"
  "$exe"
done
