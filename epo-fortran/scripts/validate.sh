#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

mode="${1:-debug}"
case "$mode" in
  debug)
    build_dir="build-validation"
    flags=(-std=f2018 -O0 -g -fcheck=all -fbacktrace -Wall -Wextra \
      -Wconversion-extra -Wimplicit-interface -Werror)
    ;;
  optimized)
    build_dir="build-validation-opt"
    flags=(-std=f2018 -O2 -Wall -Wextra -Wconversion-extra \
      -Wimplicit-interface -Werror)
    ;;
  *)
    echo "usage: scripts/validate.sh [debug|optimized]" >&2
    exit 2
    ;;
esac

rm -rf "$build_dir"
mkdir -p "$build_dir/mod" "$build_dir/obj" "$build_dir/bin"
flags+=("-J$build_dir/mod" "-I$build_dir/mod")

sources=(
  src/epo_kinds.f90
  src/epo_types.f90
  src/epo_linalg.f90
  src/epo_statistics.f90
  src/epo_core.f90
  src/epo.f90
)

objects=()
for source in "${sources[@]}"; do
  object="$build_dir/obj/$(basename "${source%.f90}").o"
  gfortran "${flags[@]}" -c "$source" -o "$object"
  objects+=("$object")
done

for source in test/*.f90 app/*.f90 example/*.f90; do
  executable="$build_dir/bin/$(basename "${source%.f90}")"
  gfortran "${flags[@]}" "$source" "${objects[@]}" -o "$executable"
  "$executable"
done
