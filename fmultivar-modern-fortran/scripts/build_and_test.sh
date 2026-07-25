#!/usr/bin/env bash
set -euo pipefail
mode=${1:-debug}
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
fc=${FC:-gfortran}
case "$mode" in
  debug)
    flags=(-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace)
    out=build/debug
    ;;
  release)
    flags=(-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace)
    out=build/release
    ;;
  *) echo "Unknown mode: $mode" >&2; exit 2 ;;
esac
mods=(
  fmultivar_kinds
  fmultivar_rng
  fmultivar_special
  fmultivar_linalg
  fmultivar_integration
  fmultivar_distributions
  fmultivar_optimizer
  fmultivar_skew
  fmultivar_grid
  fmultivar_compat
  fmultivar
)
rm -rf "$out"
mkdir -p "$out"
for mod in "${mods[@]}"; do
  "$fc" "${flags[@]}" -c "src/$mod.f90" -I "$out" -J "$out" -o "$out/$mod.o"
done
objects=()
for mod in "${mods[@]}"; do objects+=("$out/$mod.o"); done
link_flags=(-Wl,--no-warn-execstack -llapack -lblas)
for test_name in test_distributions test_utilities test_fitting test_compatibility; do
  "$fc" "${flags[@]}" "test/$test_name.f90" "${objects[@]}" -I "$out" -J "$out" \
    "${link_flags[@]}" -o "$out/$test_name"
done
for app_name in demo_fmultivar fit_csv; do
  "$fc" "${flags[@]}" "app/$app_name.f90" "${objects[@]}" -I "$out" -J "$out" \
    "${link_flags[@]}" -o "$out/$app_name"
done
"$fc" "${flags[@]}" example/integration_example.f90 "${objects[@]}" -I "$out" -J "$out" \
  "${link_flags[@]}" -o "$out/integration_example"
./test/check_license.sh
"$out/test_distributions"
"$out/test_utilities"
"$out/test_fitting"
"$out/test_compatibility"
"$out/demo_fmultivar"
"$out/integration_example"
"$out/fit_csv" data/sample_returns.csv normal
"$out/fit_csv" data/sample_returns.csv snorm 0 900
"$out/fit_csv" data/sample_returns.csv st 6 900
"$out/fit_csv" data/sample_returns.csv cauchy 0 900
echo "$mode build and test workflow passed."
