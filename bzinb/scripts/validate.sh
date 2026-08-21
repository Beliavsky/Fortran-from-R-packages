#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
rm -rf build-validation
mkdir build-validation
FF="-std=f2008 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -J build-validation -I build-validation"
SRC="src/bzinb_kinds.f90 src/bzinb_special.f90 src/bzinb_rng.f90 src/bzinb_linalg.f90 src/bzinb_optimize.f90 src/bzinb_distributions.f90 src/bzinb_em.f90 src/bzinb_fit.f90 src/bzinb.f90"
gfortran $FF -c $SRC
mv ./*.o build-validation/
gfortran $FF -c test/test_support.f90 -o build-validation/test_support.o
OBJS=$(find build-validation -maxdepth 1 -name '*.o' ! -name test_support.o | sort | tr '\n' ' ')
for src in test/test_distributions.f90 test/test_simulation.f90 test/test_poisson_fits.f90 test/test_bnb_bzinb.f90 test/test_em_parity.f90 test/test_pairwise_full.f90 test/test_weighted_pairwise.f90; do
  exe="build-validation/$(basename "${src%.f90}")"
  gfortran $FF "$src" build-validation/test_support.o $OBJS -o "$exe"
  "$exe"
done
gfortran $FF example/demo_bzinb.f90 $OBJS -o build-validation/demo_bzinb
build-validation/demo_bzinb
