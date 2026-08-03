#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
rm -rf build-gfortran
mkdir -p build-gfortran/mod build-gfortran/bin
flags="-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace"
sources="src/fixedincome_kinds.f90 src/fixedincome_types.f90 src/fixedincome_terms.f90 src/fixedincome_compounding.f90 src/fixedincome_interpolation.f90 src/fixedincome_curves.f90 src/fixedincome.f90"
for source in test/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -J build-gfortran/mod -I build-gfortran/mod $sources "$source" -o "build-gfortran/bin/$name"
  "build-gfortran/bin/$name"
done
for source in app/*.f90 example/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -J build-gfortran/mod -I build-gfortran/mod $sources "$source" -o "build-gfortran/bin/$name"
  "build-gfortran/bin/$name" >/dev/null
done
printf '%s\n' 'strict GNU Fortran validation: PASS'
