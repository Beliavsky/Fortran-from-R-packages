#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
test -f "$root/LICENSE"
grep -q "GNU GENERAL PUBLIC LICENSE" "$root/LICENSE"
grep -q 'license = "GPL-2.0-or-later"' "$root/fpm.toml"
while IFS= read -r -d '' file; do
  grep -q 'SPDX-License-Identifier: GPL-2.0-or-later' "$file"
  grep -q 'GPL-2.0-or-later\|GNU General Public License' "$file"
done < <(find "$root/src" "$root/app" "$root/example" "$root/test" -name '*.f90' -print0)
echo "GPL-2.0-or-later source license checks passed."
