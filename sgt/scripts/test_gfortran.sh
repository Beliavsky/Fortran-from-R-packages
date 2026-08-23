#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
rm -rf build-gfortran
mkdir build-gfortran
src="src/sgt_kinds.f90 src/sgt_special.f90 src/sgt_distribution.f90 src/sgt_mle.f90 src/sgt.f90"
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all -J build-gfortran -I build-gfortran -c $src
for test_src in test/*.f90; do
  name=$(basename "$test_src" .f90)
  gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all -J build-gfortran -I build-gfortran \
    ./*.o "$test_src" -o "build-gfortran/$name"
  "build-gfortran/$name"
done
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all -J build-gfortran -I build-gfortran \
  ./*.o example/basic_sgt.f90 -o build-gfortran/basic_sgt
"build-gfortran/basic_sgt"
rm -f ./*.o
