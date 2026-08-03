#!/usr/bin/env sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
OSQP_SRC="$ROOT/upstream/osqp-r/src/osqp_sources"
QDLDL_SRC="$ROOT/upstream/osqp-r/src/qdldl_sources"
OSQP_BUILD="$ROOT/backend/osqp-build"
OSQP_INSTALL="$ROOT/backend/osqp-install"
BRIDGE_BUILD="$ROOT/backend/bridge-build"
OUT_DIR="$ROOT/backend/lib"

if ! command -v cmake >/dev/null 2>&1; then
  echo "cmake is required to build the bundled OSQP backend." >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
cp "$ROOT/bridge/generic_printing.h" "$OSQP_SRC/include/private/printing.h"

cmake -S "$OSQP_SRC" -B "$OSQP_BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DOSQP_VERSION=1.0.0 \
  -DOSQP_BUILD_SHARED_LIB=OFF \
  -DOSQP_BUILD_STATIC_LIB=ON \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DOSQP_ALGEBRA_BACKEND=builtin \
  -DOSQP_ENABLE_PRINTING=ON \
  -DOSQP_ENABLE_PROFILING=ON \
  -DOSQP_ENABLE_INTERRUPT=ON \
  -DOSQP_USE_LONG=OFF \
  -DOSQP_USE_FLOAT=OFF \
  -DOSQP_BUILD_DEMO_EXE=OFF \
  -DOSQP_BUILD_UNITTESTS=OFF \
  -DOSQP_CODEGEN=OFF \
  -DOSQP_ENABLE_DERIVATIVES=OFF \
  -DFETCHCONTENT_SOURCE_DIR_QDLDL="$QDLDL_SRC" \
  -DCMAKE_INSTALL_PREFIX="$OSQP_INSTALL"
cmake --build "$OSQP_BUILD" --target install --parallel

cmake -S "$ROOT/bridge" -B "$BRIDGE_BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DOSQP_ROOT="$OSQP_INSTALL" \
  -DCMAKE_LIBRARY_OUTPUT_DIRECTORY="$OUT_DIR"
cmake --build "$BRIDGE_BUILD" --parallel

LIB=$(find "$BRIDGE_BUILD" "$OUT_DIR" -type f \( -name 'libosqp_fortran_bridge.so' -o -name 'libosqp_fortran_bridge.dylib' \) | head -n 1)
if [ -z "$LIB" ]; then
  echo "The OSQP bridge library was not produced." >&2
  exit 1
fi
DEST="$OUT_DIR/$(basename "$LIB")"
if [ "$LIB" != "$DEST" ]; then
  cp "$LIB" "$DEST"
fi

echo "Backend built: $DEST"
cd "$ROOT"
if command -v fpm >/dev/null 2>&1; then
  echo "Verifying a real QP solve through the Fortran interface..."
  fpm run --example basic_qp
fi
