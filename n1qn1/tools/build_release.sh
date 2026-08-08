#!/usr/bin/env sh
set -eu
rm -rf build-release
mkdir build-release
gfortran -std=f2018 -O3 -DNDEBUG -Wimplicit-interface \
  -Werror=implicit-interface -J build-release -I build-release \
  src/n1qn1.f90 test/test_n1qn1.f90 -o build-release/test_n1qn1
./build-release/test_n1qn1
gfortran -std=f2018 -O3 -DNDEBUG -Wimplicit-interface \
  -Werror=implicit-interface -J build-release -I build-release \
  src/n1qn1.f90 example/banana.f90 -o build-release/banana
./build-release/banana
