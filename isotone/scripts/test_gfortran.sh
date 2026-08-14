#!/usr/bin/env bash
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
build="${TMPDIR:-/tmp}/isotone-fortran-test.$$"
trap 'rm -rf "$build"' EXIT
mkdir -p "$build"
cd "$build"
flags='-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -J . -I .'
s="$root/src"
for f in \
  "$s/isotone_kinds.f90" "$s/isotone_linalg.f90" "$s/isotone_utils.f90" \
  "$s/isotone_gpava.f90" "$s/vendor_nnls/nnls_kinds.f90" \
  "$s/vendor_nnls/nnls_linalg.f90" "$s/vendor_nnls/nnls_solver.f90" \
  "$s/vendor_nnls/nnls.f90" "$s/isotone_mregnn.f90" \
  "$s/isotone_active.f90" "$s/isotone.f90"; do
  gfortran $flags -c "$f"
done
for t in "$root"/test/*.f90; do
  exe=$(basename "${t%.f90}")
  gfortran $flags "$t" ./*.o -o "$exe"
  "./$exe"
done
