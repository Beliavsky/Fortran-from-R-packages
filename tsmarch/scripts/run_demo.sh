#!/usr/bin/env sh
FLAGS=${FLAGS:-"-std=f2018 -Wall -Wextra -Werror -pedantic -O2"}
export FLAGS
. "$(dirname "$0")/common.sh"
cd "$ROOT"
compile_library
"$FC" $FLAGS -J"$MOD" -I"$MOD" "$ROOT/example/demo_tsmarch.f90" "$OBJ"/*.o -o "$BIN/demo_tsmarch"
"$BIN/demo_tsmarch"
