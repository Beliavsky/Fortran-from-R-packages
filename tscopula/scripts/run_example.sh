#!/usr/bin/env bash
set -euo pipefail
mode="${1:-opt}"
root="$(cd "$(dirname "$0")/.." && pwd)"
build="$root/build/$mode"
mod="$build/mod"
if [[ ! -f "$build/tscopula.o" ]]; then "$root/scripts/run_tests.sh" "$mode"; fi
flags=(-std=f2018 -ffree-line-length-none -Wall -Wextra -Werror -Wno-maybe-uninitialized -J"$mod" -I"$mod")
if [[ "$mode" == "check" ]]; then flags+=( -O0 -g -fcheck=all -fbacktrace ); else flags+=( -O3 ); fi
objects=()
for name in tscopula_kinds tscopula_status tscopula_math tscopula_margins tscopula_vtransforms tscopula_paircopula tscopula_timeseries tscopula_dvine tscopula_models tscopula_compat tscopula; do objects+=("$build/$name.o"); done
gfortran "${flags[@]}" "$root/example/tscopula_demo.f90" "${objects[@]}" -o "$build/tscopula_demo"
"$build/tscopula_demo"
