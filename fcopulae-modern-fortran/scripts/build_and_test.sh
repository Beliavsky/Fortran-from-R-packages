#!/usr/bin/env bash
set -euo pipefail

mode="${1:-debug}"
fc="${FC:-gfortran}"
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
rm -rf build
mkdir -p build/mod build/obj build/bin

common=(-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace)
case "$mode" in
  debug) flags=(-O0 -g -fcheck=all) ;;
  release) flags=(-O2) ;;
  *) echo "Unknown mode: $mode" >&2; exit 2 ;;
esac
link=(-llapack -lblas -Wl,--no-warn-execstack)

sources=(
  src/fcopulae_kinds.f90
  src/fcopulae_rng.f90
  src/fcopulae_special.f90
  src/fcopulae_integration.f90
  src/fcopulae_optimizer.f90
  src/fcopulae_linalg.f90
  src/fcopulae_distributions.f90
  src/fcopulae_utils.f90
  src/fcopulae_archimedean.f90
  src/fcopulae_elliptical.f90
  src/fcopulae_extreme_value.f90
  src/fcopulae_empirical.f90
  src/fcopulae.f90
)
objects=()
for source in "${sources[@]}"; do
  object="build/obj/$(basename "${source%.f90}").o"
  "$fc" "${common[@]}" "${flags[@]}" -Jbuild/mod -Ibuild/mod -c "$source" -o "$object"
  objects+=("$object")
done

"$fc" "${common[@]}" "${flags[@]}" -Jbuild/mod -Ibuild/mod -c test/test_support.f90 -o build/obj/test_support.o
for name in test_archimedean test_elliptical test_extreme_empirical; do
  "$fc" "${common[@]}" "${flags[@]}" -Jbuild/mod -Ibuild/mod "test/$name.f90" \
    build/obj/test_support.o "${objects[@]}" "${link[@]}" -o "build/bin/$name"
  "build/bin/$name"
done

for path in app/demo_fcopulae.f90 app/fit_csv.f90 example/dependence_example.f90; do
  name="$(basename "${path%.f90}")"
  "$fc" "${common[@]}" "${flags[@]}" -Jbuild/mod -Ibuild/mod "$path" "${objects[@]}" "${link[@]}" -o "build/bin/$name"
done

build/bin/demo_fcopulae
build/bin/dependence_example
build/bin/fit_csv data/sample_uv.csv archm 1 80
build/bin/fit_csv data/sample_uv.csv elliptical norm 80
build/bin/fit_csv data/sample_uv.csv ev gumbel 80
./test/check_license.sh
printf '%s build, tests, and applications passed.\n' "$mode"
