#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
fc=${FC:-gfortran}
build=build/optimized
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
flags=(-std=f2018 -Wall -Wextra -Werror -Wno-maybe-uninitialized -O3)
sources=(corpcor_kinds corpcor_types corpcor_linalg corpcor_weighted \
  corpcor_matrix_tools corpcor_shrinkage corpcor)
objects=()
for name in "${sources[@]}"; do
  "$fc" "${flags[@]}" -J "$build/mod" -I "$build/mod" -c "src/$name.f90" -o "$build/obj/$name.o"
  objects+=("$build/obj/$name.o")
done
for source in test/*.f90 example/*.f90 app/*.f90; do
  name=$(basename "${source%.f90}")
  "$fc" "${flags[@]}" -I "$build/mod" "$source" "${objects[@]}" -o "$build/bin/$name"
done
for exe in "$build"/bin/test_*; do "$exe"; done
for exe in "$build"/bin/example_*; do "$exe" >/dev/null; done
"$build/bin/demo_corpcor" >/dev/null
echo "optimized build: PASS"
