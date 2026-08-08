#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
rm -rf build-strict
mkdir build-strict
FF="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0"
SRC="src/mla_kinds.f90 src/mla_interfaces.f90 src/mla_linalg.f90 src/mla_derivatives.f90 src/mla_lmm.f90 src/marqlevalg.f90 src/marqlevalg_lmm.f90"
for t in test/*.f90; do
    exe="build-strict/$(basename "${t%.f90}")"
    gfortran $FF -J build-strict -I build-strict $SRC "$t" -o "$exe"
    "$exe"
done
for e in example/*.f90; do
    exe="build-strict/$(basename "${e%.f90}")"
    gfortran $FF -J build-strict -I build-strict $SRC "$e" -o "$exe"
    "$exe"
done
