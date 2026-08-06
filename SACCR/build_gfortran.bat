@echo off
setlocal
if "%FC%"=="" set FC=gfortran
if "%FFLAGS%"=="" set FFLAGS=-std=f2018 -O2 -Wall -Wextra

if exist build rmdir /s /q build
mkdir build

%FC% %FFLAGS% -J build -I build -c ^
 dependencies\trading\src\trading_kinds.f90 ^
 dependencies\trading\src\trading_strings.f90 ^
 dependencies\trading\src\trading_stats.f90 ^
 dependencies\trading\src\trading_dependence.f90 ^
 dependencies\trading\src\trading_dynamic_beta.f90 ^
 dependencies\trading\src\trading_trades.f90 ^
 dependencies\trading\src\trading_csa.f90 ^
 dependencies\trading\src\trading_curve.f90 ^
 dependencies\trading\src\trading_hash_table.f90 ^
 dependencies\trading\src\trading_climate.f90 ^
 dependencies\trading\src\trading_lottery.f90 ^
 dependencies\trading\src\trading_betting.f90 ^
 dependencies\trading\src\trading_io.f90 ^
 dependencies\trading\src\trading.f90 ^
 src\saccr_types.f90 ^
 src\saccr_supervisory.f90 ^
 src\saccr_core.f90 ^
 src\saccr_addon.f90 ^
 src\saccr_portfolio.f90 ^
 src\saccr_io.f90 ^
 src\saccr_examples.f90 ^
 src\saccr.f90
if errorlevel 1 exit /b 1

move /y *.o build\ >nul
%FC% %FFLAGS% -I build test\test_saccr.f90 build\*.o -o build\test_saccr.exe
if errorlevel 1 exit /b 1
build\test_saccr.exe
if errorlevel 1 exit /b 1

%FC% %FFLAGS% -I build app\saccr_demo.f90 build\*.o -o build\saccr_demo.exe
if errorlevel 1 exit /b 1
build\saccr_demo.exe
