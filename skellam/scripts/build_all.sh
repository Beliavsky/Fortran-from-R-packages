#!/usr/bin/env bash
set -euo pipefail

mode=${1:-check}
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"

case "$mode" in
  check)
    flags=(-std=f2018 -O0 -Wall -Wextra -Werror -fcheck=all -fbacktrace)
    ;;
  release)
    flags=(-std=f2018 -O3 -Wall -Wextra -Werror)
    ;;
  *)
    echo "usage: $0 [check|release]" >&2
    exit 2
    ;;
esac

build="build/$mode"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"

sources=(
  skellam_kinds
  skellam_special
  skellam_distribution
  skellam_optimization
  skellam_estimation
  skellam
)

objects=()
for source in "${sources[@]}"; do
  object="$build/obj/$source.o"
  gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" \
    -c "src/$source.f90" -o "$object"
  objects+=("$object")
done

for source in test/*.f90 example/*.f90 app/*.f90; do
  name=$(basename "$source" .f90)
  executable="$build/bin/$name"
  gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" \
    "${objects[@]}" "$source" -o "$executable"
  "$executable"
done
