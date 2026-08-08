#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build-gfortran-debug"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj"
fc=${FC:-gfortran}
flags="-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -J$build/mod -I$build/mod"
compile() {
  src=$1
  obj="$build/obj/$(basename "${src%.f90}").o"
  $fc $flags -c "$src" -o "$obj"
}
compile "$root/dependencies/lbfgs/src/lbfgs_kinds.f90"
compile "$root/dependencies/lbfgs/src/lbfgs_status.f90"
compile "$root/dependencies/lbfgs/src/lbfgs_solver.f90"
compile "$root/dependencies/lbfgs/src/lbfgs.f90"
compile "$root/dependencies/lbfgsb3/src/lbfgsb3_core.f90"
compile "$root/dependencies/lbfgsb3/src/lbfgsb3.f90"
compile "$root/src/rcppnumerical_kinds.f90"
compile "$root/src/rcppnumerical_callbacks.f90"
compile "$root/src/rcppnumerical_gk_data.f90"
compile "$root/src/rcppnumerical_integration_1d.f90"
compile "$root/src/rcppnumerical_cuhre.f90"
compile "$root/src/rcppnumerical_optimization.f90"
compile "$root/src/rcppnumerical_logistic.f90"
compile "$root/src/rcppnumerical.f90"
objects=$(find "$build/obj" -name '*.o' -print)
$fc $flags "$root/test/test_rcppnumerical.f90" $objects -o "$build/test_rcppnumerical"
"$build/test_rcppnumerical"
for example in "$root"/example/*.f90; do
  exe="$build/$(basename "${example%.f90}")"
  $fc $flags "$example" $objects -o "$exe"
  "$exe"
done
