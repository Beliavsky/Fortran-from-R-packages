#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
build="$root/build-validation"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
flags=(-std=f2018 -pedantic -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace -O0)
mv="$root/vendor/mvtnorm-fortran"
mv_sources=(
  mvtnorm_kinds mvtnorm_types mvtnorm_special mvtnorm_linalg mvtnorm_random
  mvtnorm_distributions mvtnorm_probabilities mvtnorm_triangular
  mvtnorm_conditioning mvtnorm_quantiles mvtnorm_likelihood mvtnorm
)
for src in "${mv_sources[@]}"; do
  gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" -c "$mv/src/$src.f90" -o "$build/obj/$src.o"
done
for src in matrix_normal_utils matrix_normal_distribution matrixNormal; do
  gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" -c "$root/src/$src.f90" -o "$build/obj/$src.o"
done
gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" -c "$root/test/test_support.f90" -o "$build/obj/test_support.o"
objects=("$build"/obj/mvtnorm_*.o "$build"/obj/matrix_normal_*.o "$build/obj/matrixNormal.o" "$build/obj/test_support.o")
for name in test_utils test_density test_probability test_random; do
  gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" "$root/test/$name.f90" "${objects[@]}" -o "$build/bin/$name"
  "$build/bin/$name"
done
lib_objects=("$build"/obj/mvtnorm_*.o "$build"/obj/matrix_normal_*.o "$build/obj/matrixNormal.o")
gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" "$root/example/basic.f90" "${lib_objects[@]}" -o "$build/bin/basic"
"$build/bin/basic"
