#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")" && pwd)
build="$root/build-validation"
rm -rf "$build" && mkdir -p "$build/mod" "$build/obj" "$build/bin"
flags=(-std=f2018 -pedantic -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -fbacktrace -O0)
mods=(mvtnorm_kinds mvtnorm_types mvtnorm_special mvtnorm_linalg mvtnorm_random mvtnorm_distributions mvtnorm_probabilities mvtnorm_triangular mvtnorm_conditioning mvtnorm_quantiles mvtnorm_likelihood mvtnorm)
for m in "${mods[@]}"; do gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" -c "$root/vendor/mvtnorm-fortran/src/$m.f90" -o "$build/obj/$m.o"; done
for f in "$root"/src/*.f90; do n=$(basename "$f" .f90); gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" -c "$f" -o "$build/obj/$n.o"; done
objects=("$build"/obj/*.o)
for t in "$root"/test/*.f90; do n=$(basename "$t" .f90); gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" "$t" "${objects[@]}" -o "$build/bin/$n"; "$build/bin/$n"; done
gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" "$root/example/basic.f90" "${objects[@]}" -o "$build/bin/basic"
"$build/bin/basic"
