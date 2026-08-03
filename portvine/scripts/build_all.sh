#!/usr/bin/env bash
set -euo pipefail
mode="${1:-checked}"
case "$mode" in
  checked)
    flags=(-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace)
    build=build_checked
    ;;
  optimized)
    flags=(-std=f2018 -O3 -Wall -Wextra -Wimplicit-interface)
    build=build_optimized
    ;;
  *)
    echo "usage: $0 [checked|optimized]" >&2
    exit 2
    ;;
esac
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
rugarch=(rugarch_kinds rugarch_math rugarch_rng rugarch_gh rugarch_distributions rugarch_optimizer rugarch_linalg rugarch_resampling rugarch_types rugarch_models rugarch_fit rugarch_risk rugarch_arfima rugarch_backtests rugarch_inference rugarch_evaluation rugarch_workflows rugarch_complete rugarch)
rvine=(rvine_kinds rvine_math rvine_bicop rvine_fit rvine_dvine rvine_cvine rvine_tools rvinecopulib)
portvine=(portvine_kinds portvine_types portvine_stats portvine_ordering portvine_marginals portvine_conditional portvine_dependence portvine_workflow portvine)
for m in "${rugarch[@]}"; do
  gfortran "${flags[@]}" -J "$build/mod" -I "$build/mod" -c "vendor/rugarch/src/$m.f90" -o "$build/obj/$m.o"
done
for m in "${rvine[@]}"; do
  gfortran "${flags[@]}" -J "$build/mod" -I "$build/mod" -c "vendor/rvinecopulib/src/$m.f90" -o "$build/obj/$m.o"
done
for m in "${portvine[@]}"; do
  gfortran "${flags[@]}" -J "$build/mod" -I "$build/mod" -c "src/$m.f90" -o "$build/obj/$m.o"
done
objects=("$build"/obj/*.o)
for src in test/*.f90 example/*.f90 app/*.f90; do
  exe="$build/bin/$(basename "${src%.f90}")"
  gfortran "${flags[@]}" -J "$build/mod" -I "$build/mod" "$src" "${objects[@]}" -o "$exe"
done
for exe in "$build"/bin/test_*; do "$exe"; done
for exe in "$build"/bin/example_* "$build"/bin/demo_*; do "$exe" >/dev/null; done
echo "$mode build: PASS"
