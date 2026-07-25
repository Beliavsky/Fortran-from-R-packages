#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
[ -f "$root/LICENSE" ]
grep -q 'GNU GENERAL PUBLIC LICENSE' "$root/LICENSE"
grep -q 'GPL-2.0-or-later' "$root/fpm.toml"
find "$root/src" "$root/app" "$root/example" "$root/test" -type f -name '*.f90' -print | while IFS= read -r f; do
  grep -q 'SPDX-License-Identifier: GPL-2.0-or-later' "$f" || { echo "missing SPDX header: $f"; exit 1; }
  grep -q 'under GPL version 2 or later' "$f" || { echo "missing GPL notice: $f"; exit 1; }
done
echo 'GPL-2.0-or-later source license checks passed.'
