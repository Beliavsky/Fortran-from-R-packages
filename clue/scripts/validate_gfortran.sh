#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
build=${1:-"$root/.validation-build"}
rm -rf "$build"
mkdir -p "$build"
cd "$build"
fc=${FC:-gfortran}
flags=${FFLAGS:-"-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -Werror -Wimplicit-interface -Werror=implicit-interface"}

$fc $flags -J. -I. -c \
  "$root/vendor/lpsolve-fortran/src/lpsolve_types.f90" \
  "$root/vendor/lpsolve-fortran/src/lpsolve_simplex.f90" \
  "$root/vendor/lpsolve-fortran/src/lpsolve_core.f90" \
  "$root/vendor/lpsolve-fortran/src/lpsolve_special.f90" \
  "$root/vendor/lpsolve-fortran/src/lpsolve.f90" \
  "$root/src/clue_kinds.f90" \
  "$root/src/clue_lsap.f90" \
  "$root/src/clue_partition.f90" \
  "$root/src/clue_agreement.f90" \
  "$root/src/clue_dissimilarity.f90" \
  "$root/src/clue_trees.f90" \
  "$root/src/clue_pava.f90" \
  "$root/src/clue_target_fit.f90" \
  "$root/src/clue_medoid.f90" \
  "$root/src/clue_consensus.f90" \
  "$root/src/clue_validity.f90" \
  "$root/src/clue_pclust.f90" \
  "$root/src/clue_sumt.f90" \
  "$root/src/clue_fit.f90" \
  "$root/src/clue.f90"

objs=(lpsolve_types.o lpsolve_simplex.o lpsolve_core.o lpsolve_special.o lpsolve.o \
  clue_kinds.o clue_lsap.o clue_partition.o clue_agreement.o clue_dissimilarity.o \
  clue_trees.o clue_pava.o clue_target_fit.o clue_medoid.o clue_consensus.o \
  clue_validity.o clue_pclust.o clue_sumt.o clue_fit.o clue.o)

for source in "$root"/test/*.f90; do
  name=$(basename "$source" .f90)
  $fc $flags -I. "$source" "${objs[@]}" -o "$name"
  ./$name
done
