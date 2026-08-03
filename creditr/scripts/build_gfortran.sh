#!/usr/bin/env bash
set -euo pipefail

mode="${1:-debug}"
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

case "$mode" in
  debug)
    flags=(-std=f2018 -Wall -Wextra -Werror -pedantic -fcheck=all \
      -ffpe-trap=invalid,zero,overflow -fbacktrace)
    ;;
  release)
    flags=(-std=f2018 -O3 -Wall -Wextra -Werror -pedantic)
    ;;
  *)
    echo "usage: $0 [debug|release]" >&2
    exit 2
    ;;
esac

build="build/$mode"
rm -rf "$build"
mkdir -p "$build/obj" "$build/mod" "$build/bin"

sources=(
  src/creditr_kinds.f90
  src/creditr_dates.f90
  src/creditr_curve.f90
  src/creditr_cds.f90
  src/creditr_api.f90
)

objects=()
for source in "${sources[@]}"; do
  object="$build/obj/$(basename "${source%.f90}").o"
  gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" -c "$source" -o "$object"
  objects+=("$object")
done

for source in test/*.f90 app/*.f90 example/*.f90; do
  name="$(basename "${source%.f90}")"
  exe="$build/bin/$name"
  gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" "$source" "${objects[@]}" -o "$exe"
  "$exe"
done
