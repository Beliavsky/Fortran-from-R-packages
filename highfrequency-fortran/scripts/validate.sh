#!/usr/bin/env bash
set -euo pipefail

fc="${FC:-gfortran}"
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

sources=(
  src/highfrequency_kinds.f90
  src/highfrequency_types.f90
  src/highfrequency_stats.f90
  src/highfrequency_linalg.f90
  src/highfrequency_data.f90
  src/highfrequency_cleaning.f90
  src/highfrequency_realized.f90
  src/highfrequency_optimize.f90
  src/highfrequency_models.f90
  src/highfrequency_jumps.f90
  src/highfrequency_leadlag.f90
  src/highfrequency_spot.f90
  src/highfrequency_remedi.f90
  src/highfrequency.f90
)

run_build() {
  local mode="$1"
  shift
  local flags=("$@")
  local build="build-${mode}"
  rm -rf "$build"
  mkdir -p "$build/mod" "$build/obj" "$build/bin"

  local objects=()
  for source in "${sources[@]}"; do
    local object="$build/obj/$(basename "${source%.f90}").o"
    "$fc" "${flags[@]}" -J "$build/mod" -I "$build/mod" -c "$source" -o "$object"
    objects+=("$object")
  done

  for source in test/*.f90; do
    local exe="$build/bin/$(basename "${source%.f90}")"
    "$fc" "${flags[@]}" -J "$build/mod" -I "$build/mod" "$source" "${objects[@]}" -o "$exe"
    "$exe"
  done

  for source in app/*.f90 example/*.f90; do
    local exe="$build/bin/$(basename "${source%.f90}")"
    "$fc" "${flags[@]}" -J "$build/mod" -I "$build/mod" "$source" "${objects[@]}" -o "$exe"
    "$exe" >/dev/null
  done
}

common=(-std=f2018 -Wall -Wextra -Wconversion-extra -Werror -Wno-compare-reals)
run_build debug "${common[@]}" -O0 -fcheck=all -fbacktrace
run_build optimized "${common[@]}" -O2

python3 - <<'PY'
from pathlib import Path
import hashlib
import tomllib

tomllib.loads(Path("fpm.toml").read_text())
for path in list(Path("src").glob("*.f90")) + list(Path("test").glob("*.f90")) + \
            list(Path("app").glob("*.f90")) + list(Path("example").glob("*.f90")):
    text = path.read_text()
    assert "SPDX-License-Identifier: GPL-2.0-or-later" in text, path
    assert "implicit none" in text.lower(), path
    assert all(ord(ch) < 128 for ch in text), path
    for number, line in enumerate(text.splitlines(), 1):
        assert len(line) <= 132, (path, number, len(line))
for manifest in ["provenance/original-files.sha256", "provenance/translated-files.sha256"]:
    for line in Path(manifest).read_text().splitlines():
        digest, rel = line.split("  ", 1)
        actual = hashlib.sha256(Path(rel).read_bytes()).hexdigest()
        assert actual == digest, rel
print("source_audit: PASS")
PY

echo "validation: PASS"
