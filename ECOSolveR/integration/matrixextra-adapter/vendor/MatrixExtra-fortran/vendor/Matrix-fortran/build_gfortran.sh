#!/usr/bin/env sh
set -eu

flags="-std=f2018 -O2 -Wall -Wextra -pedantic"
build=build_gfortran
rm -rf "$build"
mkdir -p "$build/mod" "$build/bin"

sources="
src/matrix_kinds.f90
src/matrix_status.f90
src/matrix_dense.f90
src/matrix_decompositions.f90
src/matrix_functions.f90
src/matrix_sparse.f90
src/matrix_sparse_solvers.f90
src/matrix_ordering.f90
src/matrix_io.f90
src/matrix_constructors.f90
src/matrix_sparse_stats.f90
src/matrix_advanced.f90
src/matrix.f90
"

objects=""
for source in $sources; do
    object="$build/$(basename "$source" .f90).o"
    gfortran $flags -J "$build/mod" -I "$build/mod" -c "$source" -o "$object"
    objects="$objects $object"
done

for source in test/*.f90 app/*.f90 example/*.f90; do
    name=$(basename "$source" .f90)
    gfortran $flags -I "$build/mod" "$source" $objects -o "$build/bin/$name"
done

for test_program in "$build"/bin/test_*; do
    "$test_program"
done

printf '%s\n' "Build complete. Programs are in $build/bin."
