#!/usr/bin/env sh
set -eu

FC="${FC:-gfortran}"
FFLAGS="${FFLAGS:--std=f2018 -O2 -Wall -Wextra}"

rm -rf build
mkdir build

trading_sources="
dependencies/trading/src/trading_kinds.f90
dependencies/trading/src/trading_strings.f90
dependencies/trading/src/trading_stats.f90
dependencies/trading/src/trading_dependence.f90
dependencies/trading/src/trading_dynamic_beta.f90
dependencies/trading/src/trading_trades.f90
dependencies/trading/src/trading_csa.f90
dependencies/trading/src/trading_curve.f90
dependencies/trading/src/trading_hash_table.f90
dependencies/trading/src/trading_climate.f90
dependencies/trading/src/trading_lottery.f90
dependencies/trading/src/trading_betting.f90
dependencies/trading/src/trading_io.f90
dependencies/trading/src/trading.f90
"

saccr_sources="
src/saccr_types.f90
src/saccr_supervisory.f90
src/saccr_core.f90
src/saccr_addon.f90
src/saccr_portfolio.f90
src/saccr_io.f90
src/saccr_examples.f90
src/saccr.f90
"

$FC $FFLAGS -J build -I build -c $trading_sources $saccr_sources
mv ./*.o build/

$FC $FFLAGS -I build test/test_saccr.f90 build/*.o -o build/test_saccr
./build/test_saccr

$FC $FFLAGS -I build app/saccr_demo.f90 build/*.o -o build/saccr_demo
./build/saccr_demo
