#!/usr/bin/env sh
set -eu
FC=${FC:-gfortran}
FLAGS=${FFLAGS:--std=f2018 -O2 -Wall -Wextra -Wpedantic}
mkdir -p build/mod build/obj build/bin
sources="src/ycevo_kinds.f90 src/ycevo_status.f90 src/ycevo_kernel.f90 src/ycevo_curve.f90 src/ycevo_types.f90 src/ycevo_linalg.f90 src/ycevo_estimation.f90 src/ycevo_prediction.f90 src/ycevo_simulation.f90 src/ycevo_io.f90 src/ycevo.f90"
objects=""
for src in $sources; do
  obj="build/obj/$(basename "${src%.f90}").o"
  $FC $FLAGS -J build/mod -I build/mod -c "$src" -o "$obj"
  objects="$objects $obj"
done
$FC $FLAGS -J build/mod -I build/mod $objects test/test_ycevo.f90 -o build/bin/test_ycevo
$FC $FLAGS -J build/mod -I build/mod $objects app/main.f90 -o build/bin/ycevo_example
