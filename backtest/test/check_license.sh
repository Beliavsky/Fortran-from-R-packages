#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)

if ! grep -q 'license = "GPL-2.0-or-later"' "$root/fpm.toml"; then
  echo "fpm.toml does not declare GPL-2.0-or-later" >&2
  exit 1
fi
if ! grep -q 'GNU GENERAL PUBLIC LICENSE' "$root/LICENSE"; then
  echo "LICENSE does not contain the GNU GPL text" >&2
  exit 1
fi
while IFS= read -r file; do
  grep -q 'SPDX-License-Identifier: GPL-2.0-or-later' "$file" || {
    echo "missing SPDX header: $file" >&2
    exit 1
  }
  grep -q 'either version 2 of the License, or any later version' "$file" || {
    echo "missing GPL-2.0-or-later notice: $file" >&2
    exit 1
  }
done < <(find "$root/src" "$root/app" "$root/example" "$root/test" -name '*.f90' -type f | sort)
echo "GPL-2.0-or-later source license checks passed."
