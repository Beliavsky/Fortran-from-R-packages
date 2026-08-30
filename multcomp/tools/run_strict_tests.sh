#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build="${TMPDIR:-/tmp}/multcomp-fortran-strict-$$"
trap 'rm -rf "$build"' EXIT
mkdir -p "$build"

flags=(
  -std=f2018
  -Wall
  -Wextra
  -Werror
  -Wimplicit-interface
  -fimplicit-none
  -fcheck=all
  -ffpe-trap=invalid,zero,overflow
  -J"$build"
  -I"$build"
)

python3 "$root/tools/check_source_rules.py"
cd "$build"

dep="$root/dependencies/mvtnorm-fortran/src"
for unit in \
  mvtnorm_kinds mvtnorm_types mvtnorm_random mvtnorm_linalg mvtnorm_special \
  mvtnorm_distributions mvtnorm_probabilities mvtnorm_triangular \
  mvtnorm_conditioning mvtnorm_quantiles mvtnorm_likelihood mvtnorm; do
  gfortran "${flags[@]}" -c "$dep/$unit.f90"
done

src="$root/src"
for unit in \
  multcomp_kinds multcomp_types multcomp_math multcomp_parm multcomp_contrasts \
  multcomp_maxsets multcomp_glht multcomp_cld multcomp_helpers multcomp_mmm multcomp; do
  gfortran "${flags[@]}" -c "$src/$unit.f90"
done

objects=("$build"/*.o)
for source in "$root"/test/*.f90; do
  name="$(basename "$source" .f90)"
  gfortran "${flags[@]}" "$source" "${objects[@]}" -o "$build/$name"
  "$build/$name"
done

for source in "$root"/example/*.f90; do
  name="$(basename "$source" .f90)"
  gfortran "${flags[@]}" "$source" "${objects[@]}" -o "$build/$name"
  "$build/$name" >/dev/null
done

echo "Strict multcomp-fortran validation passed."
