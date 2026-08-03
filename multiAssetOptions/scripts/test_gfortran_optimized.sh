#!/usr/bin/env sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build/gfortran-optimized"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"

flags="-std=f2018 -Wall -Wextra -Werror -pedantic -O3"
sources="
$root/src/mao_kinds.f90
$root/src/mao_status.f90
$root/src/mao_types.f90
$root/src/mao_sparse.f90
$root/src/mao_grid.f90
$root/src/mao_payoff.f90
$root/src/mao_operator.f90
$root/src/mao_pricing.f90
$root/src/multi_asset_options.f90
"

objects=""
for source in $sources; do
    object="$build/obj/$(basename "${source%.f90}").o"
    gfortran $flags -J "$build/mod" -I "$build/mod" -c "$source" -o "$object"
    objects="$objects $object"
done

for source in "$root"/test/*.f90 "$root"/example/*.f90 "$root"/app/*.f90; do
    name=$(basename "${source%.f90}")
    gfortran $flags -J "$build/mod" -I "$build/mod" "$source" $objects \
        -o "$build/bin/$name"
    "$build/bin/$name"
done
