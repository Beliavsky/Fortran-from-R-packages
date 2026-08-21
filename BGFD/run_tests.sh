#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
rm -rf build
mkdir -p build/mod build/obj build/bin
flags=(-std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all)
src=(
  src/adequacy_kinds.f90
  src/adequacy_interfaces.f90
  src/adequacy_math.f90
  src/adequacy_optim.f90
  src/adequacy_gof.f90
  src/adequacy_stats.f90
  src/adequacy_model.f90
  src/bgfd_core.f90
  src/bgfd_distributions.f90
  src/bgfd_fit.f90
  src/bgfd.f90
)
gfortran "${flags[@]}" -J build/mod -I build/mod -c "${src[@]}"
mv ./*.o build/obj/
for test_src in test/test_*.f90; do
  name="$(basename "${test_src%.f90}")"
  gfortran "${flags[@]}" -J build/mod -I build/mod "$test_src" build/obj/*.o -o "build/bin/$name"
  "build/bin/$name"
done
gfortran "${flags[@]}" -J build/mod -I build/mod example/example_bgfd.f90 build/obj/*.o -o build/bin/example_bgfd
build/bin/example_bgfd
