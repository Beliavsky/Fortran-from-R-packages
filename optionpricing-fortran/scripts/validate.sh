#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

sources=(
  "$root/src/optionpricing_kinds.f90"
  "$root/src/optionpricing_types.f90"
  "$root/src/optionpricing_math.f90"
  "$root/src/optionpricing_random.f90"
  "$root/src/optionpricing_linalg.f90"
  "$root/src/optionpricing_european.f90"
  "$root/src/optionpricing_asian_analytic.f90"
  "$root/src/optionpricing_asian_mc.f90"
  "$root/src/optionpricing_asian_qmc.f90"
  "$root/src/optionpricing.f90"
)
common=(-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror)

build_mode() {
  local name=$1
  shift
  local flags=("$@")
  local dir="$work/$name"
  mkdir -p "$dir/mod" "$dir/obj" "$dir/bin"

  local objects=()
  local source object
  for source in "${sources[@]}"; do
    object="$dir/obj/$(basename "${source%.f90}").o"
    gfortran "${common[@]}" "${flags[@]}" -J "$dir/mod" -I "$dir/mod" \
      -c "$source" -o "$object"
    objects+=("$object")
  done

  local target exe
  for target in "$root"/test/*.f90 "$root"/app/*.f90 "$root"/example/*.f90; do
    exe="$dir/bin/$(basename "${target%.f90}")"
    gfortran "${common[@]}" "${flags[@]}" -J "$dir/mod" -I "$dir/mod" \
      "$target" "${objects[@]}" -o "$exe"
    "$exe"
  done
  echo "$name: PASS"
}

build_mode debug -O0 -g -fcheck=all -fbacktrace
build_mode optimized -O2
python3 "$root/scripts/reference_values.py"
python3 "$root/scripts/audit.py"
echo "validation: PASS"
