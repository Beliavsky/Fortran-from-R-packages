#!/usr/bin/env sh
set -eu
rm -rf build
mkdir build
flags='-std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all'
for s in adequacy_kinds adequacy_interfaces adequacy_math adequacy_stats adequacy_optim adequacy_gof adequacy_model
do
    gfortran $flags -Jbuild -Ibuild -c "src/$s.f90" -o "build/$s.o"
done
objs='build/adequacy_kinds.o build/adequacy_interfaces.o build/adequacy_math.o build/adequacy_stats.o build/adequacy_optim.o build/adequacy_gof.o build/adequacy_model.o'
for f in test/*.f90
do
    name=$(basename "$f" .f90)
    gfortran $flags -Ibuild "$f" $objs -o "build/$name"
    "build/$name"
done
gfortran $flags -Ibuild example/example_normal.f90 $objs -o build/example_normal
build/example_normal
