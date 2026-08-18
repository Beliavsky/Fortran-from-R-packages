#!/usr/bin/env sh
set -eu
fpm build
fpm test
fpm run --example ewens_demo
