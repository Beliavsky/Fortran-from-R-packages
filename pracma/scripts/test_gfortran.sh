#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build-gfortran"
rm -rf "$build"
mkdir -p "$build"
cd "$build"
flags="-std=f2018 -Wall -Wextra -Wpedantic -Werror -Wno-compare-reals -Wno-maybe-uninitialized -Wno-intrinsic-shadow -fcheck=all -ffpe-trap=invalid,zero,overflow -O0 -g"
mods="pracma_kinds pracma_status pracma_callbacks pracma_types quadprog_kinds quadprog_core quadprog pracma_basic pracma_linalg pracma_polynomial pracma_special pracma_differentiation pracma_integration pracma_roots pracma_optimization pracma_interpolation pracma_ode pracma_signal_stats pracma_geometry pracma_combinatorics pracma_compat pracma"
for f in $mods; do
    gfortran $flags -c -J . -I . "$root/src/$f.f90"
done
objs=""
for f in $mods; do objs="$objs $f.o"; done
for source in "$root"/test/*.f90; do
    exe=$(basename "$source" .f90)
    gfortran $flags -J . -I . "$source" $objs -o "$exe"
    "./$exe"
done
for source in "$root"/example/*.f90 "$root"/app/*.f90; do
    exe=$(basename "$source" .f90)
    gfortran $flags -J . -I . "$source" $objs -o "$exe"
    "./$exe"
done
