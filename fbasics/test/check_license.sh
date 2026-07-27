#!/usr/bin/env sh
set -eu
failed=0
for f in $(find src app example test -type f -name '*.f90' | sort); do
  if ! grep -q 'SPDX-License-Identifier: GPL-2.0-or-later' "$f"; then
    echo "missing GPL-2.0-or-later SPDX header: $f" >&2
    failed=1
  fi
  if ! grep -q 'fBasics package' "$f"; then
    echo "missing fBasics attribution: $f" >&2
    failed=1
  fi
done
grep -q 'license = "GPL-2.0-or-later"' fpm.toml || failed=1
test -s LICENSE || failed=1
if [ "$failed" -ne 0 ]; then exit 1; fi
echo 'GPL-2.0-or-later source license checks passed.'
