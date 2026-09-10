#!/bin/sh
set -eu
python3 scripts/audit.py
fpm build
fpm test
fpm run --all
fpm run --example --all
