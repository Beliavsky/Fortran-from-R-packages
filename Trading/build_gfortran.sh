#!/usr/bin/env sh
set -eu

rm -rf build
mkdir build

sources="
src/trading_kinds.f90
src/trading_strings.f90
src/trading_stats.f90
src/trading_dependence.f90
src/trading_dynamic_beta.f90
src/trading_trades.f90
src/trading_csa.f90
src/trading_curve.f90
src/trading_hash_table.f90
src/trading_climate.f90
src/trading_lottery.f90
src/trading_betting.f90
src/trading_io.f90
src/trading.f90
"

gfortran -std=f2018 -O2 -J build -I build -c $sources
mv ./*.o build/
gfortran -std=f2018 -O2 -I build test/test_trading.f90 build/*.o \
  -o build/test_trading
./build/test_trading
gfortran -std=f2018 -O2 -I build app/trading_demo.f90 build/*.o \
  -o build/trading_demo
./build/trading_demo
