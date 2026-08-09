#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")/.."
rm -rf build_strict
mkdir build_strict
FLAGS="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0 -J build_strict -I build_strict"
gfortran $FLAGS -c src/dykstra_kinds.f90 -o build_strict/dykstra_kinds.o
gfortran $FLAGS -c src/dykstra_linalg.f90 -o build_strict/dykstra_linalg.o
gfortran $FLAGS -c src/dykstra_solver.f90 -o build_strict/dykstra_solver.o
gfortran $FLAGS -c src/dykstra.f90 -o build_strict/dykstra.o
OBJS="build_strict/dykstra_kinds.o build_strict/dykstra_linalg.o build_strict/dykstra_solver.o build_strict/dykstra.o"
for SRC in test/*.f90; do
    EXE="build_strict/$(basename "${SRC%.f90}")"
    gfortran $FLAGS "$SRC" $OBJS -o "$EXE"
    "$EXE"
done
for SRC in example/*.f90; do
    EXE="build_strict/$(basename "${SRC%.f90}")"
    gfortran $FLAGS "$SRC" $OBJS -o "$EXE"
    "$EXE"
done
