#!/usr/bin/env sh
set -eu
FC=${FC:-gfortran}
FLAGS="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0"
rm -rf build_gfortran
mkdir -p build_gfortran/mod build_gfortran/obj build_gfortran/bin
$FC $FLAGS -c src/cec2013_kinds.f90 -J build_gfortran/mod -o build_gfortran/obj/cec2013_kinds.o
$FC $FLAGS -c src/cec2013_data.f90 -I build_gfortran/mod -J build_gfortran/mod -o build_gfortran/obj/cec2013_data.o
$FC $FLAGS -c src/cec2013_functions.f90 -I build_gfortran/mod -J build_gfortran/mod -o build_gfortran/obj/cec2013_functions.o
$FC $FLAGS -c src/cec2013.f90 -I build_gfortran/mod -J build_gfortran/mod -o build_gfortran/obj/cec2013.o
OBJS="build_gfortran/obj/cec2013_kinds.o build_gfortran/obj/cec2013_data.o build_gfortran/obj/cec2013_functions.o build_gfortran/obj/cec2013.o"
for src in test/*.f90; do
  name=$(basename "$src" .f90)
  $FC $FLAGS "$src" $OBJS -I build_gfortran/mod -J build_gfortran/mod -o "build_gfortran/bin/$name"
  "build_gfortran/bin/$name"
done
for src in example/*.f90; do
  name=$(basename "$src" .f90)
  $FC $FLAGS "$src" $OBJS -I build_gfortran/mod -J build_gfortran/mod -o "build_gfortran/bin/$name"
  "build_gfortran/bin/$name" >/dev/null
done
