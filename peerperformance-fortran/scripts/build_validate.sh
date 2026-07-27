#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
build=${1:-build-validation}
flags=${FFLAGS:--O0 -g -std=f2018 -Wall -Wextra -Werror -Wconversion-extra -fcheck=all -fbacktrace}
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
modules=(
  peerperformance_kinds peerperformance_types peerperformance_math
  peerperformance_linalg peerperformance_stats peerperformance_pi
  peerperformance_bootstrap peerperformance_screening peerperformance
)
objects=()
for name in "${modules[@]}"; do
  gfortran $flags -J"$build/mod" -I"$build/mod" -c "src/$name.f90" -o "$build/obj/$name.o"
  objects+=("$build/obj/$name.o")
done
for source in test/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -I"$build/mod" "$source" "${objects[@]}" -o "$build/bin/$name"
  "$build/bin/$name"
done
for source in app/*.f90 example/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -I"$build/mod" "$source" "${objects[@]}" -o "$build/bin/$name"
  "$build/bin/$name" >/dev/null
done
echo 'validation: PASS'
