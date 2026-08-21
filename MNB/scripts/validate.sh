#!/usr/bin/env sh
set -eu
rm -rf build objects
mkdir -p build objects
FFLAGS="-std=f2008 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -O0"
SRC="src/mnb_kinds.f90 src/mnb_types.f90 src/mnb_math.f90 src/mnb_optimizer.f90 src/mnb_core.f90 src/mnb_simulation.f90 src/mnb_residuals.f90 src/mnb_influence.f90 src/mnb_envelope.f90 src/mnb.f90"
gfortran $FFLAGS -J build -I build -c $SRC
mv ./*.o objects/
for t in test/*.f90; do
  exe="build/$(basename "$t" .f90)"
  gfortran $FFLAGS -J build -I build objects/*.o "$t" -o "$exe"
  "$exe"
done
gfortran $FFLAGS -J build -I build objects/*.o example/demo_mnb.f90 -o build/demo_mnb
build/demo_mnb
