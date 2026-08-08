#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
build="$root/build-gfortran-opt"
rm -rf "$build"
mkdir -p "$build"
cd "$build"

flags=(-std=f2018 -O3 -march=native -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface)
sources=(
  "$root/dependencies/lbfgsb3/src/lbfgsb3_core.f90"
  "$root/dependencies/lbfgsb3/src/lbfgsb3.f90"
  "$root/src/nfcp_types.f90"
  "$root/src/nfcp_math.f90"
  "$root/src/nfcp_parameters.f90"
  "$root/src/nfcp_kalman.f90"
  "$root/src/nfcp_forecast.f90"
  "$root/src/nfcp_simulation.f90"
  "$root/src/nfcp_analysis.f90"
  "$root/src/nfcp_stitch.f90"
  "$root/src/nfcp_options.f90"
  "$root/src/nfcp_mle.f90"
  "$root/src/nfcp.f90"
)

gfortran "${flags[@]}" -J . -I . -c "${sources[@]}"
gfortran "${flags[@]}" -I . "$root/test/test_nfcp.f90" ./*.o -o test_nfcp
./test_nfcp
gfortran "${flags[@]}" -I . "$root/example/two_factor_oil.f90" ./*.o -o two_factor_oil
./two_factor_oil
