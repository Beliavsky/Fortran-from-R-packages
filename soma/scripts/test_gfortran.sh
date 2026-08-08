#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
rm -rf build_strict
mkdir build_strict
cd build_strict
flags=(-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all)
gfortran "${flags[@]}" -J . -I . -c \
  ../src/soma_kinds.f90 \
  ../src/soma_types.f90 \
  ../src/soma_random.f90 \
  ../src/soma_optimizer.f90 \
  ../src/soma.f90
objs=(soma_kinds.o soma_types.o soma_random.o soma_optimizer.o soma.o)
for source in ../test/*.f90; do
  exe="$(basename "${source%.f90}")"
  gfortran "${flags[@]}" -J . -I . "${objs[@]}" "$source" -o "$exe"
  "./$exe"
done
for source in ../example/*.f90; do
  exe="$(basename "${source%.f90}")"
  gfortran "${flags[@]}" -J . -I . "${objs[@]}" "$source" -o "$exe"
  "./$exe"
done
