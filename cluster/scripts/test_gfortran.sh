#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build-gfortran"
rm -rf "$build"
mkdir -p "$build"
cd "$build"
flags="-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow"
src="$root/src/fastcluster_kinds.f90 $root/src/fastcluster_types.f90 $root/src/fastcluster_distances.f90 $root/src/fastcluster_core.f90 $root/src/cluster_types.f90 $root/src/cluster_linalg.f90 $root/src/cluster_daisy.f90 $root/src/cluster_partition.f90 $root/src/cluster_hierarchy.f90 $root/src/cluster_diagnostics.f90 $root/src/cluster_ellipsoid.f90 $root/src/cluster.f90"
# shellcheck disable=SC2086
gfortran $flags -J . -c $src
ar rcs libcluster.a ./*.o
for file in "$root"/test/*.f90; do
  exe=$(basename "$file" .f90)
  # shellcheck disable=SC2086
  gfortran $flags -I . "$file" libcluster.a -o "$exe"
  "./$exe"
done
for file in "$root"/example/*.f90 "$root"/app/*.f90; do
  exe=$(basename "$file" .f90)
  # shellcheck disable=SC2086
  gfortran $flags -I . "$file" libcluster.a -o "$exe"
  "./$exe" >/dev/null
done
