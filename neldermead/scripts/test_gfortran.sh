#!/usr/bin/env bash
set -euo pipefail
FC=${FC:-gfortran}
FLAGS="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0"
rm -rf build-strict
mkdir -p build-strict/mod build-strict/obj build-strict/bin
sources=(src/neldermead_kinds.f90 src/neldermead_types.f90 src/neldermead_simplex.f90 src/neldermead_core.f90 src/neldermead_frontends.f90 src/neldermead.f90)
objs=()
for s in "${sources[@]}"; do
  o="build-strict/obj/$(basename "${s%.f90}").o"
  $FC $FLAGS -Jbuild-strict/mod -Ibuild-strict/mod -c "$s" -o "$o"
  objs+=("$o")
done
for t in test/*.f90; do
  exe="build-strict/bin/$(basename "${t%.f90}")"
  $FC $FLAGS -Jbuild-strict/mod -Ibuild-strict/mod "$t" "${objs[@]}" -o "$exe"
  "$exe"
done
for e in example/*.f90; do
  exe="build-strict/bin/$(basename "${e%.f90}")"
  $FC $FLAGS -Jbuild-strict/mod -Ibuild-strict/mod "$e" "${objs[@]}" -o "$exe"
  "$exe"
done
