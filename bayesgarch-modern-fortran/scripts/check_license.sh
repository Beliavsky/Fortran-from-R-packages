#!/bin/sh
# SPDX-License-Identifier: GPL-2.0-or-later
set -eu

fail=0
for file in $(find src app example test -type f -name '*.f90' | sort); do
    if ! grep -q 'SPDX-License-Identifier: GPL-2.0-or-later' "$file"; then
        echo "missing GPL-2.0-or-later SPDX identifier: $file" >&2
        fail=1
    fi
    if ! grep -q 'Copyright (C) 2008-2021 David Ardia' "$file"; then
        echo "missing original copyright notice: $file" >&2
        fail=1
    fi
    if ! grep -q 'either version 2 of the License' "$file"; then
        echo "missing GPL version 2-or-later notice: $file" >&2
        fail=1
    fi
done

if ! grep -q '^license = "GPL-2.0-or-later"' fpm.toml; then
    echo "fpm.toml does not declare GPL-2.0-or-later" >&2
    fail=1
fi
if ! grep -q 'GNU GENERAL PUBLIC LICENSE' LICENSE; then
    echo "LICENSE does not contain the GNU GPL text" >&2
    fail=1
fi
if ! grep -q 'Version 2, June 1991' LICENSE; then
    echo "LICENSE is not the GPL version 2 text" >&2
    fail=1
fi

if [ "$fail" -ne 0 ]; then
    exit 1
fi

echo "GPL-2.0-or-later source license checks passed."
