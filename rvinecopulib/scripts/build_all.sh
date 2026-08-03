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
mods=(rvine_kinds rvine_math rvine_bicop rvine_fit rvine_dvine rvine_cvine rvine_tools rvinecopulib)
for m in "${mods[@]}"; do
  gfortran "${flags[@]}" -J "$build/mod" -I "$build/mod" -c "src/$m.f90" -o "$build/obj/$m.o"
done
objects=("$build"/obj/*.o)
for src in test/*.f90 example/*.f90 app/*.f90; do
  exe="$build/bin/$(basename "${src%.f90}")"
  gfortran "${flags[@]}" -J "$build/mod" -I "$build/mod" "$src" "${objects[@]}" -o "$exe"
done
for exe in "$build"/bin/test_*; do
  "$exe"
done
for exe in "$build"/bin/example_* "$build"/bin/demo_*; do
  "$exe" >/dev/null
done
echo "$mode build: PASS"
