#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
build="$root/build-gfortran"
rm -rf "$build"
mkdir -p "$build"
cd "$build"
flags=(-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0)
for f in \
  "$root/vendor/splines-fortran-v0.1.0/src/splines_kinds.f90" \
  "$root/vendor/splines-fortran-v0.1.0/src/splines_linalg.f90" \
  "$root/vendor/splines-fortran-v0.1.0/src/splines_core.f90" \
  "$root/vendor/splines-fortran-v0.1.0/src/splines_basis.f90" \
  "$root/vendor/splines-fortran-v0.1.0/src/splines.f90" \
  "$root/src/survival_kinds.f90" "$root/src/survival_types.f90" \
  "$root/src/survival_linalg.f90" "$root/src/survival_nonparametric.f90" \
  "$root/src/survival_cox.f90" "$root/src/survival_aft.f90" \
  "$root/src/survival_stats.f90" "$root/src/survival_utils.f90" \
  "$root/src/survival_pspline.f90" "$root/src/survival.f90"; do
  gfortran "${flags[@]}" -c "$f"
done
objs=( *.o )
for testsrc in "$root"/test/*.f90; do
  exe="$(basename "${testsrc%.f90}")"
  gfortran "${flags[@]}" "$testsrc" "${objs[@]}" -o "$exe"
  "./$exe"
done
for exsrc in "$root"/example/*.f90; do
  exe="$(basename "${exsrc%.f90}")"
  gfortran "${flags[@]}" "$exsrc" "${objs[@]}" -o "$exe"
  "./$exe"
done
