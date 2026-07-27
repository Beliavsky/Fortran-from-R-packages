#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
mapfile -t files < <(find "$root/src" "$root/app" "$root/example" "$root/test" -type f -name '*.f90' | sort)
if [[ ${#files[@]} -eq 0 ]]; then
  echo "No Fortran files found" >&2
  exit 1
fi
for file in "${files[@]}"; do
  grep -q 'SPDX-License-Identifier: GPL-2.0-or-later' "$file" || {
    echo "Missing SPDX header: $file" >&2
    exit 1
  }
  grep -q 'either version 2 of the License, or' "$file" || {
    echo "Missing GPL-2-or-later notice: $file" >&2
    exit 1
  }
  grep -q 'Original garchx package copyright' "$file" || {
    echo "Missing original-package attribution: $file" >&2
    exit 1
  }
done
grep -q 'license = "GPL-2.0-or-later"' "$root/fpm.toml"
grep -q 'GNU GENERAL PUBLIC LICENSE' "$root/LICENSE"
grep -q 'Version 2, June 1991' "$root/LICENSE"
echo "GPL-2.0-or-later source license checks passed."
