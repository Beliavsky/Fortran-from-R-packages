#!/usr/bin/env sh
# SPDX-License-Identifier: GPL-3.0-only
set -eu

make clean
make check
make clean
make all fit_csv
./build/demo_rmgarch
./build/fit_csv data/sample_returns.csv
