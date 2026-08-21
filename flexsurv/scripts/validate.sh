#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
B="${TMPDIR:-/tmp}/flexsurv_fortran_validate"
rm -rf "$B"; mkdir -p "$B"; cd "$B"
FLAGS=(-std=f2008 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -J. -I.)
SRC=(flexsurv_kinds.f90 flexsurv_math.f90 flexsurv_distributions.f90 flexsurv_splines2ns.f90 flexsurv_spline.f90 \
quadprog_kinds.f90 quadprog_core.f90 quadprog.f90 \
numderiv_kinds.f90 numderiv_types.f90 numderiv_callbacks.f90 numderiv_core.f90 numderiv.f90 \
survival_kinds.f90 survival_types.f90 survival_linalg.f90 survival_aft.f90 \
desolve_kinds.f90 desolve_types.f90 desolve_rk.f90 relsurv_kinds.f90 relsurv_ratetable.f90 \
flexsurv_fit.f90 flexsurv_spline_fit.f90 flexsurv_spline_interactions.f90 \
flexsurv_mixture.f90 flexsurv_mixture_full.f90 flexsurv_fmixmsm.f90 \
flexsurv_multistate.f90 flexsurv_multistate_uncertainty.f90 flexsurv_shared_multistate.f90 flexsurv_final_states.f90 \
flexsurv_standardize.f90 flexsurv_standardize_advanced.f90 flexsurv_ajfit.f90 \
flexsurv_rtrunc.f90 flexsurv_fracpoly.f90 flexsurv_diagnostics.f90 flexsurv_custom.f90 flexsurv.f90)
for f in "${SRC[@]}"; do gfortran "${FLAGS[@]}" -c "$ROOT/src/$f"; done
OBJS=( ./*.o )
for t in "$ROOT"/test/*.f90; do
  exe="$(basename "$t" .f90)"
  gfortran "${FLAGS[@]}" "$t" "${OBJS[@]}" -o "$exe"
  "./$exe"
done
gfortran "${FLAGS[@]}" "$ROOT/example/demo_flexsurv.f90" "${OBJS[@]}" -o demo_flexsurv
./demo_flexsurv
