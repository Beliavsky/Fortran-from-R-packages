#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build_gfortran"
rm -rf "$build"
mkdir -p "$build"
cd "$root"
flags="-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -J$build -I$build"
sources="src/fmbasics_kinds.f90 src/fmbasics_dates.f90 src/fmbasics_interpolation.f90 src/fmbasics_conventions.f90 src/fmbasics_rates.f90 src/fmbasics_curves.f90 src/fmbasics_credit.f90 src/fmbasics_volatility.f90 src/fmbasics_money.f90 src/fmbasics.f90"
objects=""
for source in $sources; do
  object="$build/$(basename "$source" .f90).o"
  gfortran $flags -c "$source" -o "$object"
  objects="$objects $object"
done
for source in test/*.f90; do
  executable="$build/$(basename "$source" .f90)"
  gfortran $flags "$source" $objects -o "$executable"
  "$executable"
done
for source in example/*.f90 app/*.f90; do
  executable="$build/$(basename "$source" .f90)"
  gfortran $flags "$source" $objects -o "$executable"
  "$executable"
done
