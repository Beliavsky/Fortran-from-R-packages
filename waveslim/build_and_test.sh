#!/usr/bin/env sh
set -eu
make clean
make check
make release
make app example
./build/check/app/waveslim_demo
./build/check/example/image_packet_example
