#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build-gfortran"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
flags="-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow"
sources="greeks_kinds greeks_types greeks_math greeks_payoffs greeks_integrals greeks_black_scholes greeks_binomial greeks_monte_carlo greeks_api greeks"
objects=""
for name in $sources; do
  gfortran $flags -J"$build/mod" -I"$build/mod" -c "$root/src/$name.f90" -o "$build/obj/$name.o"
  objects="$objects $build/obj/$name.o"
done
for source in "$root"/test/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -J"$build/mod" -I"$build/mod" "$source" $objects -o "$build/bin/$name"
  "$build/bin/$name"
done
