#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
set -euo pipefail
mode="${1:-debug}"
make MODE="$mode" clean
make MODE="$mode" check
