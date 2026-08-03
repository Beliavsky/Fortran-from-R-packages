#!/usr/bin/env sh
set -eu

FC=${FC:-gfortran}
FLAGS=${FFLAGS:-"-std=f2018 -Wall -Wextra -Werror -pedantic -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -g"}
BUILD=build-gfortran
rm -rf "$BUILD"
mkdir "$BUILD"
cd "$BUILD"

for source in \
  ../src/icsnp_kinds.f90 \
  ../src/icsnp_status.f90 \
  ../src/icsnp_types.f90 \
  ../src/icsnp_linalg.f90 \
  ../src/icsnp_special.f90 \
  ../src/icsnp_pairs.f90 \
  ../src/icsnp_estimators.f90 \
  ../src/icsnp_tests.f90 \
  ../src/icsnp.f90
do
  $FC $FLAGS -J. -I. -c "$source"
done

OBJECTS="icsnp_kinds.o icsnp_status.o icsnp_types.o icsnp_linalg.o icsnp_special.o icsnp_pairs.o icsnp_estimators.o icsnp_tests.o icsnp.o"

for source in ../test/*.f90
do
  name=$(basename "$source" .f90)
  $FC $FLAGS -J. -I. $OBJECTS "$source" -o "$name"
  ./"$name"
done

for source in ../example/*.f90 ../app/*.f90
do
  name=$(basename "$source" .f90)
  $FC $FLAGS -J. -I. $OBJECTS "$source" -o "$name"
  ./"$name"
done
