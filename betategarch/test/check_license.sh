#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
[ -f "$root/LICENSE" ]
grep -q "GNU GENERAL PUBLIC LICENSE" "$root/LICENSE"
grep -q 'license = "GPL-2.0-only"' "$root/fpm.toml"
find "$root/src" "$root/app" "$root/test" "$root/example" -type f -name '*.f90' | while IFS= read -r file; do
  grep -q 'SPDX-License-Identifier: GPL-2.0-only' "$file" || { echo "Missing SPDX header: $file"; exit 1; }
  grep -q 'Copyright (C) 2011-2025 Genaro Sucarrat' "$file" || { echo "Missing original copyright: $file"; exit 1; }
  grep -q 'terms of the GNU General Public License version 2 only' "$file" || { echo "Missing GPL-2-only notice: $file"; exit 1; }
done
echo "GPL-2.0-only source license checks passed."
