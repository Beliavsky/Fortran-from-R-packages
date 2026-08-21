#!/usr/bin/env sh
set -eu
rm -rf build-direct
mkdir build-direct
rm -f ./*.o
flags="-std=f2008 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -O0 -g"
gfortran $flags -J build-direct -I build-direct -c \
  src/survival_kinds.f90 src/survival_types.f90 src/survival_linalg.f90 src/survival_cox.f90 \
  src/relsurv_kinds.f90 src/relsurv_ratetable.f90 src/relsurv_parsers.f90 \
  src/mstate_kinds.f90 src/mstate_types.f90 src/mstate_transitions.f90 src/mstate_data.f90 src/mstate_crprep.f90 \
  src/mstate_msfit.f90 src/mstate_cox.f90 src/mstate_redrank.f90 src/mstate_probtrans.f90 \
  src/mstate_nonparametric.f90 src/mstate_simulation.f90 src/mstate_relative.f90 src/mstate_utilities.f90 \
  src/mstate_markov.f90 src/mstate.f90
for t in test/*.f90; do
  exe="build-direct/$(basename "${t%.f90}")"
  gfortran $flags -J build-direct -I build-direct "$t" ./*.o -o "$exe"
  "$exe"
done
for e in example/*.f90; do
  exe="build-direct/$(basename "${e%.f90}")"
  gfortran $flags -J build-direct -I build-direct "$e" ./*.o -o "$exe"
  "$exe"
done
