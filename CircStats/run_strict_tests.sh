#!/bin/sh
set -eu
rm -rf build
mkdir build
rm -f ./*.o
srcs="src/circstats_kinds.f90 src/circstats_types.f90 src/circstats_special.f90 src/circstats_utils.f90 src/circstats_rao_table.f90 src/circstats_core.f90 src/circstats_distributions.f90 src/circstats_tests.f90 src/circstats_models.f90 src/circstats.f90"
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all -J build -I build -c $srcs
for t in test/*.f90; do
    exe="build/$(basename "$t" .f90)"
    gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all -J build -I build ./*.o "$t" -o "$exe"
    "$exe"
done
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all -J build -I build ./*.o example/basic.f90 -o build/basic
./build/basic
rm -f ./*.o
