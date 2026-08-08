#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
build="$root/build-strict"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
fc=${FC:-gfortran}
flags=(-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -J"$build/mod" -I"$build/mod")
for src in rceim_kinds rceim_random rceim_utils rceim_benchmarks rceim; do
  "$fc" "${flags[@]}" -c "$root/src/$src.f90" -o "$build/obj/$src.o"
done
objs=("$build/obj/rceim_kinds.o" "$build/obj/rceim_random.o" "$build/obj/rceim_utils.o" "$build/obj/rceim_benchmarks.o" "$build/obj/rceim.o")
for t in "$root"/test/*.f90; do
  name=$(basename "$t" .f90)
  "$fc" "${flags[@]}" "$t" "${objs[@]}" -o "$build/bin/$name"
  "$build/bin/$name"
done
for e in "$root"/example/*.f90; do
  name=$(basename "$e" .f90)
  "$fc" "${flags[@]}" "$e" "${objs[@]}" -o "$build/bin/$name"
  "$build/bin/$name"
done
