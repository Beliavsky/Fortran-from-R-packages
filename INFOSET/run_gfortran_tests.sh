#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build/gfortran-debug"
rm -rf "$build"
mkdir -p "$build"
cd "$build"
flags='-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace'
sources='quadprog_kinds quadprog_core quadprog infoset_kinds infoset_status infoset_types infoset_stats infoset_mixture infoset_core infoset_portfolio infoset'
for name in $sources; do
  gfortran $flags -J. -c "$root/src/$name.f90"
done
objects='quadprog_kinds.o quadprog_core.o quadprog.o infoset_kinds.o infoset_status.o infoset_types.o infoset_stats.o infoset_mixture.o infoset_core.o infoset_portfolio.o infoset.o'
for source in "$root"/test/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -I. "$source" $objects -o "$name"
  "./$name"
done
for source in "$root"/example/*.f90 "$root"/app/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -I. "$source" $objects -o "$name"
  "./$name"
done
