#!/usr/bin/env sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
rm -rf build-validation
mkdir build-validation

flags="-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface \
-Werror -fcheck=all -fbacktrace -O0 -g -J build-validation -I build-validation"

objects=""
for name in portfoliooptim_kinds portfoliooptim_types portfoliooptim_linalg \
  portfoliooptim_simplex portfoliooptim_risk portfoliooptim_benders \
  portfoliooptim_projection portfoliooptim
do
  gfortran $flags -c "src/$name.f90" -o "build-validation/$name.o"
  objects="$objects build-validation/$name.o"
done

for source in test/test_*.f90
do
  name=$(basename "$source" .f90)
  gfortran $flags "$source" $objects -o "build-validation/$name"
  "build-validation/$name"
done

for source in app/*.f90 example/*.f90
do
  name=$(basename "$source" .f90)
  gfortran $flags "$source" $objects -o "build-validation/$name"
  "build-validation/$name" >/dev/null
done

echo "validation: PASS"
