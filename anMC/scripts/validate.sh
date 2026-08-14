#!/usr/bin/env sh
set -eu
FC=${FC:-gfortran}
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BUILD="$ROOT/build-validation"
rm -rf "$BUILD"
mkdir -p "$BUILD"
cd "$BUILD"
FLAGS="-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -J. -I."
for f in \
  anmc_kinds.f90 anmc_types.f90 anmc_utils.f90 anmc_math.f90 \
  anmc_sampling.f90 anmc_active.f90 anmc_mc.f90 \
  anmc_probabilities.f90 anmc_conservative.f90 anmc.f90
do
  "$FC" $FLAGS -c "$ROOT/src/$f"
done
for t in test_core test_active test_orthant test_conservative
do
  "$FC" $FLAGS "$ROOT/test/$t.f90" ./*.o -o "$t"
  "./$t"
done
"$FC" $FLAGS "$ROOT/example/equicorrelated_orthant.f90" ./*.o -o example_equicorrelated_orthant
./example_equicorrelated_orthant
