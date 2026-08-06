#!/usr/bin/env sh
set -eu

FC=${FC:-gfortran}
AR=${AR:-ar}
MODE=${WQC_BUILD_MODE:-check}

if [ "$MODE" = "release" ]; then
  FLAGS="-std=f2018 -O3 -Wall -Wextra -Werror"
else
  FLAGS="-std=f2018 -O0 -g -Wall -Wextra -Werror -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow"
fi

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
BUILD="$ROOT/build/$MODE"
MOD="$BUILD/mod"
OBJ="$BUILD/obj"
BIN="$BUILD/bin"
rm -rf "$BUILD"
mkdir -p "$MOD" "$OBJ" "$BIN"

compile() {
  "$FC" $FLAGS -J "$MOD" -I "$MOD" -c "$1" -o "$OBJ/$(basename "${1%.f90}").o"
}

for name in qcsis_kinds qcsis_statistics qcsis; do
  compile "$ROOT/dependencies/qcsis/src/$name.f90"
done

for name in waveslim_kinds waveslim_status waveslim_types waveslim_math \
  waveslim_linalg waveslim_filters waveslim_transform_1d waveslim_packet \
  waveslim_transform_nd waveslim_dualtree waveslim_hilbert_stats \
  waveslim_statistics waveslim_denoise waveslim_long_memory \
  waveslim_extended waveslim; do
  compile "$ROOT/dependencies/waveslim/src/$name.f90"
done

for name in wqc_kinds wqc_random wqc_statistics wqc; do
  compile "$ROOT/src/$name.f90"
done

"$AR" rcs "$BUILD/libwqc.a" "$OBJ"/*.o
"$FC" $FLAGS -J "$MOD" -I "$MOD" "$ROOT/test/test_wqc.f90" \
  "$BUILD/libwqc.a" -o "$BIN/test_wqc"
"$BIN/test_wqc"
"$FC" $FLAGS -J "$MOD" -I "$MOD" "$ROOT/app/wqc_demo.f90" \
  "$BUILD/libwqc.a" -o "$BIN/wqc_demo"

echo "Built $BIN/wqc_demo"
