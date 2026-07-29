#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
build="$root/build-optimized"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
flags=(-std=f2018 -pedantic -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -O2)
sources=(mvtnorm_kinds mvtnorm_types mvtnorm_special mvtnorm_linalg mvtnorm_random mvtnorm_distributions mvtnorm_probabilities mvtnorm_triangular mvtnorm_conditioning mvtnorm_quantiles mvtnorm_likelihood mvtnorm)
for src in "${sources[@]}"; do
  gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" -c "$root/src/$src.f90" -o "$build/obj/$src.o"
done
gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" -c "$root/test/test_support.f90" -o "$build/obj/test_support.o"
objects=("$build"/obj/mvtnorm_*.o "$build/obj/test_support.o")
for test in "$root"/test/test_*.f90; do
  name=$(basename "$test" .f90); [[ "$name" == test_support ]] && continue
  gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" "$test" "${objects[@]}" -o "$build/bin/$name"
  "$build/bin/$name"
done
for source in "$root"/app/*.f90 "$root"/example/*.f90; do
  name=$(basename "$source" .f90)
  gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" "$source" "${objects[@]}" -o "$build/bin/$name"
  "$build/bin/$name" >/dev/null
done
