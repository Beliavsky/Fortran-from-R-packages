#!/usr/bin/env bash
set -euo pipefail
mode="${1:-checked}"
fc="${FC:-gfortran}"
case "$mode" in
  checked)
    flags="-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Werror -fcheck=all -fbacktrace"
    ;;
  optimized)
    flags="-std=f2018 -O3 -Wall -Wextra -Wpedantic -Werror"
    ;;
  *) echo "usage: $0 [checked|optimized]" >&2; exit 2 ;;
esac
root="$(cd "$(dirname "$0")/.." && pwd)"
build="$root/build/$mode"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
sources=(
  spantest_kinds spantest_types spantest_probability spantest_random
  spantest_linalg spantest_classical spantest_gl spantest_as
  spantest_simulation spantest
)
for name in "${sources[@]}"; do
  "$fc" $flags -J"$build/mod" -I"$build/mod" -c "$root/src/$name.f90" -o "$build/obj/$name.o"
done
for test_src in "$root"/test/*.f90; do
  name="$(basename "${test_src%.f90}")"
  "$fc" $flags -J"$build/mod" -I"$build/mod" "$test_src" "$build"/obj/*.o -o "$build/bin/$name"
  "$build/bin/$name"
done
"$fc" $flags -J"$build/mod" -I"$build/mod" "$root/example/spantest_demo.f90" \
  "$build"/obj/*.o -o "$build/bin/spantest_demo"
"$build/bin/spantest_demo"
