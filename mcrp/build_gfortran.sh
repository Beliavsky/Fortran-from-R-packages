#!/usr/bin/env bash
set -euo pipefail
mode="${1:-strict}"
case "$mode" in
  strict)
    flags=(-std=f2018 -O0 -g -Wall -Wextra -Werror -pedantic \
      -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    ;;
  release)
    flags=(-std=f2018 -O3 -Wall -Wextra -Werror -pedantic)
    ;;
  *)
    echo "usage: $0 [strict|release]" >&2
    exit 2
    ;;
esac
build="build/$mode"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
compile() {
  local src="$1" obj="$2"
  gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" -c "$src" -o "$obj"
}
compile src/mcrp_kinds.f90 "$build/obj/mcrp_kinds.o"
compile src/mcrp_moments.f90 "$build/obj/mcrp_moments.o"
compile src/mcrp_optimizer.f90 "$build/obj/mcrp_optimizer.o"
compile src/mcrp.f90 "$build/obj/mcrp.o"
objects=("$build/obj/mcrp_kinds.o" "$build/obj/mcrp_moments.o" \
  "$build/obj/mcrp_optimizer.o" "$build/obj/mcrp.o")
for target in test/test_moments.f90 test/test_optimizer.f90 \
  app/mcrp_demo.f90 example/moment_example.f90; do
  name="$(basename "$target" .f90)"
  gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" \
    "$target" "${objects[@]}" -o "$build/bin/$name"
done
"$build/bin/test_moments"
"$build/bin/test_optimizer"
"$build/bin/mcrp_demo" >/dev/null
"$build/bin/moment_example" >/dev/null
echo "mcrp $mode build and tests passed"
