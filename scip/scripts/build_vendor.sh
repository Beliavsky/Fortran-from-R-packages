#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
PKG="$ROOT/original/scip-master"
SCIP_SRC="$PKG/inst/scip"
SOPLEX_SRC="$PKG/inst/soplex"
CONFIG="$PKG/inst/config"
BUILD="$ROOT/vendor/build"
LIB="$ROOT/vendor/lib"
JOBS=${JOBS:-2}
BUILD_TYPE=${SCIP_FORTRAN_BUILD_TYPE:-Release}
TPI=${SCIP_FORTRAN_TPI:-none}

mkdir -p "$BUILD" "$LIB"
created_dirs=()
make_stub() {
    local d=$1
    if [ ! -d "$d" ]; then
        mkdir -p "$d"
        printf '# generated build stub\n' > "$d/CMakeLists.txt"
        created_dirs+=("$d")
    fi
}
cleanup() {
    local d
    for d in "${created_dirs[@]:-}"; do rm -rf "$d"; done
}
trap cleanup EXIT

make_stub "$SOPLEX_SRC/check"
for d in check tests doc examples applications; do make_stub "$SCIP_SRC/$d"; done

EXTRA_RELEASE=()
if [ "${SCIP_FORTRAN_FAST_BUILD:-0}" = "1" ]; then
    EXTRA_RELEASE+=("-DCMAKE_C_FLAGS_RELEASE=-O0 -DNDEBUG")
    EXTRA_RELEASE+=("-DCMAKE_CXX_FLAGS_RELEASE=-O0 -DNDEBUG")
fi

cmake -S "$SOPLEX_SRC" -B "$BUILD/soplex" \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DBUILD_SHARED_LIBS=OFF \
    -DZLIB=OFF -DGMP=OFF -DMPFR=OFF -DBOOST=OFF -DPAPILO=OFF -DQUADMATH=OFF \
    "-DCMAKE_CXX_FLAGS=-I$CONFIG" "${EXTRA_RELEASE[@]}"
cmake --build "$BUILD/soplex" --target libsoplex -j "$JOBS"

cmake -S "$SCIP_SRC" -B "$BUILD/scip" \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DBUILD_SHARED_LIBS=OFF \
    -DSOPLEX_DIR="$BUILD/soplex" \
    -DZLIB=OFF -DGMP=OFF -DREADLINE=OFF -DPAPILO=OFF -DZIMPL=OFF -DAMPL=OFF \
    -DIPOPT=OFF -DWORHP=OFF -DCONOPT=OFF -DLAPACK=OFF -DAUTOBUILD=OFF \
    -DSHARED=OFF -DEXPRINT=none -DLPS=spx -DTPI="$TPI" \
    "-DCMAKE_C_FLAGS=-I$CONFIG -DSCIP_LONGINT_FORMAT=\\\"lld\\\"" \
    "-DCMAKE_CXX_FLAGS=-I$CONFIG -DSCIP_LONGINT_FORMAT=\\\"lld\\\"" \
    "${EXTRA_RELEASE[@]}"
cmake --build "$BUILD/scip" --target libscip -j "$JOBS"

CC=${CC:-gcc}
CXX=${CXX:-g++}
AR=${AR:-ar}
"$CC" -c "$ROOT/csrc/scip_fortran_shim.c" -o "$BUILD/scip_fortran_shim.o" \
    -std=c11 -O2 -fPIC -I"$SCIP_SRC/src" -I"$BUILD/scip" -I"$CONFIG"
"$CXX" -c "$ROOT/csrc/standalone_streams.cpp" -o "$BUILD/standalone_streams.o" \
    -std=c++17 -O2 -fPIC -I"$CONFIG"

rm -f "$LIB/libscipfortran_backend.a"
if [ "$(uname -s)" = "Darwin" ]; then
    libtool -static -o "$LIB/libscipfortran_backend.a" \
        "$BUILD/scip/lib/libscip.a" "$BUILD/soplex/lib/libsoplex.a" \
        "$BUILD/scip_fortran_shim.o" "$BUILD/standalone_streams.o"
else
    cat > "$BUILD/merge.mri" <<MRI
create $LIB/libscipfortran_backend.a
addlib $BUILD/scip/lib/libscip.a
addlib $BUILD/soplex/lib/libsoplex.a
addmod $BUILD/scip_fortran_shim.o
addmod $BUILD/standalone_streams.o
save
end
MRI
    "$AR" -M < "$BUILD/merge.mri"
fi

echo "Built $LIB/libscipfortran_backend.a"
