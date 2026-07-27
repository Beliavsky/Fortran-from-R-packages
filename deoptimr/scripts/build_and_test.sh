#!/usr/bin/env bash
set -euo pipefail
mode=${1:-debug}
root=$(cd "$(dirname "$0")/.." && pwd)
fc=${FC:-gfortran}
common=(-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace)
case "$mode" in
    debug) flags=("${common[@]}" -O0 -g -fcheck=all) ;;
    release) flags=("${common[@]}" -O2) ;;
    *) echo "usage: $0 debug|release" >&2; exit 2 ;;
esac
build="$root/build/$mode"
obj="$build/obj"
mod="$build/mod"
bin="$build/bin"
rm -rf "$build"
mkdir -p "$obj" "$mod" "$bin"

sources=(
    deoptimr_kinds
    deoptimr_interfaces
    deoptimr_rng
    deoptimr_types
    deoptimr_utils
    deoptimr_jde
    deoptimr_ncde
    deoptimr
)
objects=()
for name in "${sources[@]}"; do
    "$fc" "${flags[@]}" -J"$mod" -I"$mod" -c "$root/src/$name.f90" -o "$obj/$name.o"
    objects+=("$obj/$name.o")
done
"$fc" "${flags[@]}" -J"$mod" -I"$mod" -c "$root/test/test_support.f90" -o "$obj/test_support.o"
"$fc" "${flags[@]}" -J"$mod" -I"$mod" -c "$root/test/test_benchmarks.f90" -o "$obj/test_benchmarks.o"

"$fc" "${flags[@]}" -J"$mod" -I"$mod" "$root/test/test_jde.f90" \
    "${objects[@]}" "$obj/test_support.o" "$obj/test_benchmarks.o" -o "$bin/test_jde"
"$fc" "${flags[@]}" -J"$mod" -I"$mod" "$root/test/test_ncde.f90" \
    "${objects[@]}" "$obj/test_support.o" "$obj/test_benchmarks.o" -o "$bin/test_ncde"
"$fc" "${flags[@]}" -J"$mod" -I"$mod" "$root/app/demo_deoptimr.f90" \
    "${objects[@]}" -o "$bin/demo_deoptimr"
"$fc" "${flags[@]}" -J"$mod" -I"$mod" "$root/example/constrained_example.f90" \
    "${objects[@]}" -o "$bin/constrained_example"

"$bin/test_jde"
"$bin/test_ncde"
"$bin/demo_deoptimr"
"$bin/constrained_example"
"$root/test/check_license.sh"
echo "$mode build, tests, and applications passed."
