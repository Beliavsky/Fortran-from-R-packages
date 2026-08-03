#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
mode=${1:-strict}
build="$root/build-$mode"
rm -rf "$build"
mkdir -p "$build/mod" "$build/bin"

common=(-std=f2018 -J"$build/mod" -I"$build/mod")
if [[ "$mode" == "strict" ]]; then
    flags=(-O0 -g -Wall -Wextra -Werror -pedantic -fcheck=all \
        -ffpe-trap=invalid,zero,overflow -fbacktrace)
elif [[ "$mode" == "release" ]]; then
    flags=(-O3 -Wall -Wextra -Werror -pedantic)
else
    echo "usage: $0 [strict|release]" >&2
    exit 2
fi

gfortran "${common[@]}" "${flags[@]}" -c "$root/src/corpmetrics.f90" -o "$build/corpmetrics.o"

for source in "$root"/test/*.f90 "$root"/app/*.f90 "$root"/example/*.f90; do
    name=$(basename "${source%.f90}")
    gfortran "${common[@]}" "${flags[@]}" "$build/corpmetrics.o" "$source" -o "$build/bin/$name"
done

for exe in "$build"/bin/test_*; do
    "$exe"
done

"$build/bin/demo_corpmetrics" >/dev/null
"$build/bin/investment_example" >/dev/null
"$build/bin/loan_example" >/dev/null

echo "corpmetrics $mode build: PASS"
