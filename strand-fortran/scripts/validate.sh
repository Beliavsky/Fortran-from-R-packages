#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"

modules=(
  strand_kinds strand_types strand_linalg strand_stats strand_simplex
  strand_optimizer strand_data strand_simulation strand
)

build_configuration() {
  local name=$1
  shift
  local flags=("$@")
  local dir="build/validation-$name"
  rm -rf "$dir"
  mkdir -p "$dir/mod" "$dir/obj" "$dir/bin"

  for module in "${modules[@]}"; do
    gfortran -c "${flags[@]}" -J "$dir/mod" -I "$dir/mod" \
      "src/${module}.f90" -o "$dir/obj/${module}.o"
  done

  local objects=("$dir"/obj/*.o)
  for source in test/*.f90; do
    local target
    target=$(basename "$source" .f90)
    gfortran "${flags[@]}" -J "$dir/mod" -I "$dir/mod" "$source" \
      "${objects[@]}" -o "$dir/bin/$target"
    "$dir/bin/$target"
  done

  for source in app/*.f90 example/*.f90; do
    local target
    target=$(basename "$source" .f90)
    gfortran "${flags[@]}" -J "$dir/mod" -I "$dir/mod" "$source" \
      "${objects[@]}" -o "$dir/bin/$target"
    "$dir/bin/$target"
  done
}

build_configuration checked \
  -std=f2018 -O0 -g -Wall -Wextra -Wconversion-extra \
  -Wimplicit-interface -Werror -fcheck=all -fbacktrace

build_configuration optimized \
  -std=f2018 -O2 -Wall -Wextra -Wconversion-extra \
  -Wimplicit-interface -Werror

python scripts/audit_release.py
printf '%s\n' 'validation: PASS'
