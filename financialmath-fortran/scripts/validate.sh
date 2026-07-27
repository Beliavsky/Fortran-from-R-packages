#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
rm -rf build-validation
a=(src/financialmath_kinds.f90 src/financialmath_types.f90 \
   src/financialmath_math.f90 src/financialmath_cashflows.f90 \
   src/financialmath_annuities.f90 src/financialmath_loans.f90 \
   src/financialmath_derivatives.f90 src/financialmath.f90)
for opt in "-O0 -fcheck=all -fbacktrace" "-O2"; do
   rm -f ./*.o ./*.mod
   mkdir -p build-validation
   flags="-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror $opt"
   gfortran $flags -J build-validation -I build-validation -c "${a[@]}" test/test_support.f90
   for name in test_cashflows test_annuities test_loans_bonds test_derivatives; do
      gfortran $flags -J build-validation -I build-validation test/$name.f90 ./*.o -o build-validation/$name
      build-validation/$name
   done
   for source in app/financialmath_demo.f90 example/cashflows_and_loans.f90 example/options_and_forwards.f90; do
      exe=build-validation/$(basename "$source" .f90)
      gfortran $flags -J build-validation -I build-validation "$source" ./*.o -o "$exe"
      "$exe" >/dev/null
   done
done
rm -f ./*.o ./*.mod
printf '%s\n' 'validation: PASS'
