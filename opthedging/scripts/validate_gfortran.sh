#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
build="$root/build-validation"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"

opt_level=${OPT_LEVEL:--O0}
flags=(-std=f2018 "$opt_level" -g -Wall -Wextra -Wconversion-extra -Wimplicit-interface \
  -Werror -fcheck=all -fbacktrace -ffree-line-length-132 -J"$build/mod" -I"$build/mod")

sources=(
  opthedging_kinds.f90
  opthedging_interpolation.f90
  opthedging_statistics.f90
  opthedging_types.f90
  opthedging_iid.f90
  opthedging_rng.f90
  opthedging.f90
)

objects=()
for source in "${sources[@]}"; do
  object="$build/obj/${source%.f90}.o"
  gfortran "${flags[@]}" -c "$root/src/$source" -o "$object"
  objects+=("$object")
done

for source in "$root"/test/*.f90; do
  name=$(basename "${source%.f90}")
  gfortran "${flags[@]}" "$source" "${objects[@]}" -o "$build/bin/$name"
  "$build/bin/$name"
done

for source in "$root"/app/*.f90 "$root"/example/*.f90; do
  name=$(basename "${source%.f90}")
  gfortran "${flags[@]}" "$source" "${objects[@]}" -o "$build/bin/$name"
  "$build/bin/$name" >/dev/null
done

echo "validation: PASS"
