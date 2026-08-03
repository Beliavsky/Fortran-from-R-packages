#!/usr/bin/env sh
set -eu

FC=${FC:-gfortran}
BUILD=build-gfortran
rm -rf "$BUILD"
mkdir -p "$BUILD/mod" "$BUILD/obj" "$BUILD/bin"

STRICT="-std=f2008 -pedantic -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow"
SRC="src/sgt_kinds.f90 src/sgt_status.f90 src/sgt_types.f90 src/sgt_linalg.f90 src/sgt_operators.f90 src/sgt_utils.f90 src/sgt_updates.f90 src/sgt_objectives.f90 src/sgt_initial_graph.f90 src/sgt_spectral_learning.f90 src/sgt_other_learning.f90 src/spectral_graph_topology.f90"

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

"$FC" $STRICT -J "$BUILD/mod" -I "$BUILD/mod" app/demo_spectral_graph_topology.f90 $OBJECTS -o "$BUILD/bin/demo"
"$BUILD/bin/demo"

for example_source in example/*.f90; do
  name=$(basename "${example_source%.f90}")
  "$FC" $STRICT -J "$BUILD/mod" -I "$BUILD/mod" "$example_source" $OBJECTS -o "$BUILD/bin/$name"
  "$BUILD/bin/$name"
done

# Optimized library compile and smoke test.
mkdir -p "$BUILD/opt-mod" "$BUILD/opt-obj"
OPT="-std=f2008 -O3 -Wall -Wextra -Wimplicit-interface -Werror"
for source in $SRC; do
  object="$BUILD/opt-obj/$(basename "${source%.f90}").o"
  "$FC" $OPT -J "$BUILD/opt-mod" -I "$BUILD/opt-mod" -c "$source" -o "$object"
done
OPT_OBJECTS=$(find "$BUILD/opt-obj" -name '*.o' -print | sort)
"$FC" $OPT -J "$BUILD/opt-mod" -I "$BUILD/opt-mod" app/demo_spectral_graph_topology.f90 $OPT_OBJECTS -o "$BUILD/bin/demo-optimized"
"$BUILD/bin/demo-optimized" >/dev/null
