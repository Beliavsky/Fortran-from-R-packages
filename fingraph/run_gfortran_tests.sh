#!/usr/bin/env sh
set -eu

FC=${FC:-gfortran}
BUILD=build-gfortran
rm -rf "$BUILD"
mkdir -p "$BUILD/mod" "$BUILD/obj" "$BUILD/bin"

STRICT="-std=f2008 -pedantic -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow"
SRC="src/fingraph_kinds.f90 src/fingraph_status.f90 src/fingraph_types.f90 src/fingraph_linalg.f90 src/fingraph_operators.f90 src/fingraph_utils.f90 src/fingraph_learning.f90 src/fingraph_rng.f90 src/fingraph.f90"

for source in $SRC; do
  object="$BUILD/obj/$(basename "${source%.f90}").o"
  "$FC" $STRICT -J "$BUILD/mod" -I "$BUILD/mod" -c "$source" -o "$object"
done

OBJECTS=$(find "$BUILD/obj" -name '*.o' -print | sort)
for test_source in test/*.f90; do
  name=$(basename "${test_source%.f90}")
  "$FC" $STRICT -J "$BUILD/mod" -I "$BUILD/mod" "$test_source" $OBJECTS -o "$BUILD/bin/$name"
  "$BUILD/bin/$name"
done

for app_source in app/*.f90; do
  name=$(basename "${app_source%.f90}")
  "$FC" $STRICT -J "$BUILD/mod" -I "$BUILD/mod" "$app_source" $OBJECTS -o "$BUILD/bin/$name"
  "$BUILD/bin/$name"
done

for example_source in example/*.f90; do
  name=$(basename "${example_source%.f90}")
  "$FC" $STRICT -J "$BUILD/mod" -I "$BUILD/mod" "$example_source" $OBJECTS -o "$BUILD/bin/$name"
  "$BUILD/bin/$name"
done

mkdir -p "$BUILD/opt-mod" "$BUILD/opt-obj"
OPT="-std=f2008 -O3 -Wall -Wextra -Wimplicit-interface -Werror"
for source in $SRC; do
  object="$BUILD/opt-obj/$(basename "${source%.f90}").o"
  "$FC" $OPT -J "$BUILD/opt-mod" -I "$BUILD/opt-mod" -c "$source" -o "$object"
done
OPT_OBJECTS=$(find "$BUILD/opt-obj" -name '*.o' -print | sort)
for test_source in test/*.f90; do
  name=$(basename "${test_source%.f90}")
  "$FC" $OPT -J "$BUILD/opt-mod" -I "$BUILD/opt-mod" "$test_source" $OPT_OBJECTS -o "$BUILD/bin/$name-optimized"
  "$BUILD/bin/$name-optimized"
done
