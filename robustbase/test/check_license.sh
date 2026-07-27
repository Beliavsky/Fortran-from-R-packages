#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
status=0
while IFS= read -r file; do
  grep -q 'SPDX-License-Identifier: GPL-2.0-or-later' "$file" || { echo "missing SPDX: $file"; status=1; }
  grep -q 'GPL version 2 or later' "$file" || { echo "missing GPL notice: $file"; status=1; }
done < <(find "$root/src" "$root/app" "$root/example" "$root/test" -name '*.f90' -type f | sort)
grep -q 'license = "GPL-2.0-or-later"' "$root/fpm.toml" || { echo 'wrong fpm license'; status=1; }
test -s "$root/LICENSE" || { echo 'missing LICENSE'; status=1; }
if [[ $status -ne 0 ]]; then exit $status; fi
echo 'GPL-2.0-or-later source license checks passed.'
