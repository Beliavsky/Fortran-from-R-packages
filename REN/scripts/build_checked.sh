#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
fc=${FC:-gfortran}
build=build/checked
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
flags=(-std=f2018 -Wall -Wextra -Wpedantic -Werror -Wno-maybe-uninitialized -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -O0 -g)

corpcor_sources=(corpcor_kinds corpcor_types corpcor_linalg corpcor_weighted corpcor_matrix_tools \
  corpcor_shrinkage corpcor)
ren_sources=(ren_kinds ren_types ren_linalg ren_random ren_regularization ren_portfolio ren_analysis ren)
objects=()
for name in "${corpcor_sources[@]}"; do
  "$fc" "${flags[@]}" -J "$build/mod" -I "$build/mod" -c "vendor/corpcor/src/$name.f90" -o "$build/obj/$name.o"
  objects+=("$build/obj/$name.o")
done
for name in "${ren_sources[@]}"; do
  "$fc" "${flags[@]}" -J "$build/mod" -I "$build/mod" -c "src/$name.f90" -o "$build/obj/$name.o"
  objects+=("$build/obj/$name.o")
done
ar rcs "$build/libren.a" "${objects[@]}"
for source in test/*.f90; do
  name=$(basename "${source%.f90}")
  "$fc" "${flags[@]}" -I "$build/mod" "$source" "$build/libren.a" -o "$build/bin/$name"
  "$build/bin/$name"
done
for source in example/*.f90 app/*.f90; do
  name=$(basename "${source%.f90}")
  "$fc" "${flags[@]}" -I "$build/mod" "$source" "$build/libren.a" -o "$build/bin/$name"
done
"$build/bin/basic_portfolios" >/dev/null
"$build/bin/demo_ren" >/dev/null
echo "checked GNU Fortran build: PASS"
