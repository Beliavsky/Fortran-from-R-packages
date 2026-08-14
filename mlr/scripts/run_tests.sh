#!/usr/bin/env sh
set -eu
fpm test
fpm run --example basic_workflow
