#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="${TMPDIR:-/tmp}/mclust-fortran-validation-$$"
trap 'rm -rf "$build"' EXIT INT TERM
mkdir -p "$build/mod" "$build/obj" "$build/bin"
legacy="$root/src/legacy"
for f in mclust_legacy mclust2_legacy util_legacy; do
  gfortran -O2 -c "$legacy/$f.f90" -J "$build/mod" -o "$build/obj/$f.o"
done
gfortran -std=f2018 -O2 -c "$legacy/legacy_support.f90" \
  -J "$build/mod" -o "$build/obj/legacy_support.o"
flags='-std=f2018 -O2 -fcheck=all -Wall -Wextra -Werror -Wno-maybe-uninitialized -Wimplicit-interface'
mods='mclust_kinds mclust_types mclust_legacy_interfaces mclust_math mclust_linalg mclust_utilities mclust_crimcoords mclust_hierarchical mclust_models mclust_selection mclust_simulation mclust_density mclust_impute mclust_dr mclust_combine mclust_weighted mclust_bootstrap mclust_ssc mclust_classification mclust_api'
for m in $mods; do
  gfortran $flags -I "$build/mod" -J "$build/mod" -c "$root/src/$m.f90" \
    -o "$build/obj/$m.o"
done
objs=$(find "$build/obj" -name '*.o' -print | tr '\n' ' ')
for t in "$root"/test/*.f90; do
  name=$(basename "$t" .f90)
  gfortran $flags -I "$build/mod" "$t" $objs -llapack -lblas -o "$build/bin/$name"
  "$build/bin/$name"
done
gfortran $flags -I "$build/mod" "$root/example/basic_mclust.f90" $objs \
  -llapack -lblas -o "$build/bin/basic_mclust"
"$build/bin/basic_mclust"
cd "$root"
gfortran $flags -I "$build/mod" "$root/validation/diabetes_reference.f90" $objs \
  -llapack -lblas -o "$build/bin/diabetes_reference"
"$build/bin/diabetes_reference"
