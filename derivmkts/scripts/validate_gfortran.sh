#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
build="$root/build-validation"
rm -rf "$build"
mkdir -p "$build/mod" "$build/bin"

flags=(-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface \
       -Werror -fcheck=all -fbacktrace)

sources=(
  derivmkts_kinds.f90
  derivmkts_types.f90
  derivmkts_math.f90
  derivmkts_rng.f90
  derivmkts_black_scholes.f90
  derivmkts_implied.f90
  derivmkts_bonds.f90
  derivmkts_asian_analytic.f90
  derivmkts_barriers.f90
  derivmkts_perpetual.f90
  derivmkts_compound.f90
  derivmkts_jumps.f90
  derivmkts_binomial.f90
  derivmkts_simulation.f90
  derivmkts_asian_mc.f90
  derivmkts_greeks.f90
  derivmkts_quincunx.f90
  derivmkts.f90
)

objects=()
for source in "${sources[@]}"; do
  object="$build/${source%.f90}.o"
  gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" \
    -c "$root/src/$source" -o "$object"
  objects+=("$object")
done
ar rcs "$build/libderivmkts.a" "${objects[@]}"

for source in "$root"/test/*.f90 "$root"/app/*.f90 "$root"/example/*.f90; do
  name=$(basename "${source%.f90}")
  gfortran "${flags[@]}" -I"$build/mod" "$source" \
    "$build/libderivmkts.a" -o "$build/bin/$name"
  "$build/bin/$name"
done
