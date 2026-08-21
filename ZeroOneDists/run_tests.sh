#!/usr/bin/env sh
set -eu
FC=${FC:-gfortran}
FFLAGS=${FFLAGS:--std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all}
rm -rf build
mkdir -p build/mod build/obj build/bin
for f in src/zero_one_kinds.f90 src/zero_one_special.f90 src/zero_one_distributions.f90 src/zero_one_families.f90 src/zero_one_fit.f90 src/zero_one_dists.f90; do
  obj="build/obj/$(basename "${f%.f90}").o"
  $FC $FFLAGS -Jbuild/mod -Ibuild/mod -c "$f" -o "$obj"
done
objs="build/obj/zero_one_kinds.o build/obj/zero_one_special.o build/obj/zero_one_distributions.o build/obj/zero_one_families.o build/obj/zero_one_fit.o build/obj/zero_one_dists.o"
for t in test/*.f90; do
  exe="build/bin/$(basename "${t%.f90}")"
  $FC $FFLAGS -Ibuild/mod "$t" $objs -o "$exe"
  "$exe"
done
$FC $FFLAGS -Ibuild/mod example/basic.f90 $objs -o build/bin/basic
build/bin/basic
