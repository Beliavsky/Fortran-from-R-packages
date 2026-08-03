#!/usr/bin/env bash
set -euo pipefail
mode="${1:-strict}"
root="$(cd "$(dirname "$0")/.." && pwd)"
build="$root/build-$mode"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
if [[ "$mode" == strict ]]; then
  flags=(-std=f2018 -O0 -g -Wall -Wextra -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow)
else
  flags=(-std=f2018 -O3 -Wall -Wextra -Werror)
fi
sources=(esback_kinds esback_types esback_math esback_optimizer esback_esreg esback_backtests esback)
objects=()
for src in "${sources[@]}"; do
  obj="$build/obj/$src.o"
  gfortran "${flags[@]}" -J "$build/mod" -I "$build/mod" -c "$root/src/$src.f90" -o "$obj"
  objects+=("$obj")
done
for file in "$root"/test/*.f90 "$root"/app/*.f90 "$root"/example/*.f90; do
  name="$(basename "$file" .f90)"
  gfortran "${flags[@]}" -I "$build/mod" "$file" "${objects[@]}" -llapack -lblas -o "$build/bin/$name"
done
for exe in "$build"/bin/test_*; do "$exe"; done
"$build/bin/main" >/dev/null
"$build/bin/all_backtests" >/dev/null
"$build/bin/bootstrap_esr" >/dev/null
