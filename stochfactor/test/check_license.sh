#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
test -f "$root/LICENSE"
grep -q 'GNU GENERAL PUBLIC LICENSE' "$root/LICENSE"
grep -q 'license = "GPL-2.0-or-later"' "$root/fpm.toml"
while IFS= read -r f; do
  grep -q 'SPDX-License-Identifier: GPL-2.0-or-later' "$f"
  grep -q 'either version 2, or (at your option) any later version' "$f"
done < <(find "$root/src" "$root/app" "$root/example" "$root/test" -name '*.f90' -type f)
echo 'GPL-2.0-or-later source license checks passed.'
