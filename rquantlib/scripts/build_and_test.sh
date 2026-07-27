#!/usr/bin/env bash
set -euo pipefail
mode="${1:-debug}"
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
rm -rf "build/$mode"
mkdir -p "build/$mode/mod" "build/$mode/obj" "build/$mode/bin"
case "$mode" in
  debug)
    flags=(-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace)
    ;;
  release)
    flags=(-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace)
    ;;
  *)
    echo "mode must be debug or release" >&2
    exit 2
    ;;
esac
sources=(
  src/rq_kinds.f90
  src/rq_math.f90
  src/rq_optimize.f90
  src/rq_linalg.f90
  src/rq_dates.f90
  src/rq_options.f90
  src/rq_curves.f90
  src/rq_bonds.f90
  src/rq_sabr.f90
  src/rq_hull_white.f90
)
objs=()
for src in "${sources[@]}"; do
  base="$(basename "${src%.f90}")"
  obj="build/$mode/obj/$base.o"
  gfortran "${flags[@]}" -J "build/$mode/mod" -I "build/$mode/mod" -c "$src" -o "$obj"
  objs+=("$obj")
done
gfortran "${flags[@]}" -J "build/$mode/mod" -I "build/$mode/mod" \
  -c test/test_support.f90 -o "build/$mode/obj/test_support.o"
test_obj="build/$mode/obj/test_support.o"
for test_name in test_options test_dates_curves_bonds test_sabr_hull_white; do
  gfortran "${flags[@]}" -J "build/$mode/mod" -I "build/$mode/mod" \
    "test/$test_name.f90" "${objs[@]}" "$test_obj" -llapack -lblas \
    -o "build/$mode/bin/$test_name"
  "build/$mode/bin/$test_name"
done
for app_name in demo_rquantlib price_option fit_curve price_bond; do
  gfortran "${flags[@]}" -J "build/$mode/mod" -I "build/$mode/mod" \
    "app/$app_name.f90" "${objs[@]}" -llapack -lblas \
    -o "build/$mode/bin/$app_name"
done
gfortran "${flags[@]}" -J "build/$mode/mod" -I "build/$mode/mod" \
  example/curve_and_option_example.f90 "${objs[@]}" -llapack -lblas \
  -o "build/$mode/bin/curve_and_option_example"
"build/$mode/bin/demo_rquantlib"
"build/$mode/bin/price_option" european call 100 100 0.01 0.04 1 0.2
"build/$mode/bin/price_bond" 100 0.04 5 2 0.03
"build/$mode/bin/fit_curve" data/example_yield_curve.csv ns >/dev/null
"build/$mode/bin/curve_and_option_example"
"$root/test/check_license.sh"
echo "$mode build, tests, and applications passed."
