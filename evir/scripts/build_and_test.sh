#!/bin/sh
# SPDX-License-Identifier: GPL-2.0-or-later
set -eu

mode=${1:-debug}
fc=${FC:-gfortran}
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build/$mode"
obj="$build/obj"
mod="$build/mod"
bin="$build/bin"

case "$mode" in
    debug)
        flags="-std=f2018 -Wall -Wextra -Werror -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace"
        ;;
    release)
        flags="-std=f2018 -Wall -Wextra -Werror -O3"
        ;;
    *)
        echo "usage: $0 [debug|release]" >&2
        exit 2
        ;;
esac

rm -rf "$build"
mkdir -p "$obj" "$mod" "$bin"

sources="
evir_kinds.f90
evir_types.f90
evir_math.f90
evir_optimize.f90
evir_distributions.f90
evir_data.f90
evir_fitting.f90
evir_eda.f90
evir_bivariate.f90
evir.f90
"

for source in $sources; do
    base=${source%.f90}
    "$fc" $flags -J "$mod" -I "$mod" -c "$root/src/$source" -o "$obj/$base.o"
done

objects=$(find "$obj" -name '*.o' -type f | sort | tr '\n' ' ')
"$fc" $flags -J "$mod" -I "$mod" -c "$root/test/test_support.f90" -o "$obj/test_support.o"
test_objects="$objects $obj/test_support.o"

for test in test_distributions test_eda test_fitting test_bivariate; do
    "$fc" $flags -J "$mod" -I "$mod" $test_objects "$root/test/$test.f90" -o "$bin/$test"
    "$bin/$test"
done

"$fc" $flags -J "$mod" -I "$mod" $objects "$root/app/demo_evir.f90" -o "$bin/demo_evir"
"$fc" $flags -J "$mod" -I "$mod" $objects "$root/example/block_maxima_example.f90" -o "$bin/block_maxima_example"
"$bin/demo_evir" >/dev/null
"$bin/block_maxima_example" >/dev/null

echo "evir-fortran $mode build and tests passed"
