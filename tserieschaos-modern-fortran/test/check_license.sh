#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
[[ -s "$root/LICENSE" ]]
grep -q 'GNU GENERAL PUBLIC LICENSE' "$root/LICENSE"
grep -q 'Version 2, June 1991' "$root/LICENSE"
grep -q 'license = "GPL-2.0-only"' "$root/fpm.toml"
while IFS= read -r file; do
  grep -q '^! SPDX-License-Identifier: GPL-2.0-only$' "$file"
  grep -q 'GNU General Public License version 2 only' "$file"
done < <(find "$root/src" "$root/app" "$root/example" "$root/test" -type f -name '*.f90' | sort)
echo 'GPL-2.0-only source license checks passed.'
