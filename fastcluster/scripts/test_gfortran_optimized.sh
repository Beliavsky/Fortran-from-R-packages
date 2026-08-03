#!/bin/sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build-gfortran-optimized"
rm -rf "$build"
mkdir -p "$build/mod" "$build/bin"

flags="-std=f2018 -Wall -Wextra -Wpedantic -Werror -O3"

compile_module() {
  src=$1
  obj=$2
  gfortran $flags -J "$build/mod" -I "$build/mod" -c "$root/$src" -o "$build/$obj"
}

compile_module src/fastcluster_kinds.f90 fastcluster_kinds.o
compile_module src/fastcluster_types.f90 fastcluster_types.o
compile_module src/fastcluster_distances.f90 fastcluster_distances.o
compile_module src/fastcluster_core.f90 fastcluster_core.o
compile_module src/fastcluster.f90 fastcluster.o

objects="$build/fastcluster_kinds.o $build/fastcluster_types.o $build/fastcluster_distances.o $build/fastcluster_core.o $build/fastcluster.o"

for source in "$root"/test/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -J "$build/mod" -I "$build/mod" "$source" $objects -o "$build/bin/$name"
  "$build/bin/$name"
done

for source in "$root"/app/*.f90 "$root"/example/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -J "$build/mod" -I "$build/mod" "$source" $objects -o "$build/bin/$name"
  "$build/bin/$name" >/dev/null
done
