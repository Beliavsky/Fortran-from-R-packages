#!/usr/bin/env sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
HIGHS_SRC="$ROOT/upstream/highs-r/inst/HiGHS"
HIGHS_BUILD="$ROOT/backend/highs-build"
HIGHS_INSTALL="$ROOT/backend/highs-install"
BRIDGE_BUILD="$ROOT/backend/bridge-build"
OUT_DIR="$ROOT/backend/lib"

if ! command -v cmake >/dev/null 2>&1; then
  echo "cmake is required to build the bundled HiGHS backend." >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
export CPLUS_INCLUDE_PATH="$ROOT/bridge/compat${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"

cmake -S "$HIGHS_SRC" -B "$HIGHS_BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$HIGHS_INSTALL" \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DBUILD_SHARED_LIBS=OFF \
  -DBUILD_TESTING=OFF \
  -DBUILD_EXAMPLES=OFF \
  -DZLIB=OFF \
  -DFAST_BUILD=ON \
  -DCUPDLP_GPU=OFF
cmake --build "$HIGHS_BUILD" --target install --parallel

cmake -S "$ROOT/bridge" -B "$BRIDGE_BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DHIGHS_ROOT="$HIGHS_INSTALL" \
  -DCMAKE_LIBRARY_OUTPUT_DIRECTORY="$OUT_DIR"
cmake --build "$BRIDGE_BUILD" --parallel

LIB=$(find "$BRIDGE_BUILD" "$OUT_DIR" -type f \( -name 'libhighs_fortran_bridge.so' -o -name 'libhighs_fortran_bridge.dylib' \) | head -n 1)
if [ -z "$LIB" ]; then
  echo "The bridge library was not produced." >&2
  exit 1
fi
cp "$LIB" "$OUT_DIR/$(basename "$LIB")"

echo "Backend built: $OUT_DIR/$(basename "$LIB")"
echo "Run from the package root with: fpm run"
