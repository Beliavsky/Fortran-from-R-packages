#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
build=${TMPDIR:-/tmp}/mixsqp-fortran-strict
rm -rf "$build"
mkdir -p "$build"
flags=(-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -J"$build" -I"$build")
src=(mixsqp_kinds mixsqp_lapack mixsqp_types mixsqp_utils mixsqp_em mixsqp_solver mixsqp_highlevel mixsqp_simulate mixsqp)
for base in "${src[@]}"; do
  gfortran "${flags[@]}" -c "$root/src/$base.f90" -o "$build/$base.o"
done
objs=("$build"/*.o)
for file in "$root"/test/*.f90; do
  exe="$build/$(basename "${file%.f90}")"
  gfortran "${flags[@]}" "$file" "${objs[@]}" -llapack -lblas -o "$exe"
  "$exe"
done
gfortran "${flags[@]}" "$root/example/basic.f90" "${objs[@]}" -llapack -lblas -o "$build/basic"
"$build/basic"
