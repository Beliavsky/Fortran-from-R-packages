@echo off
setlocal
if not exist build mkdir build
del /q build\*.o build\*.mod build\test_xva.exe build\xva_demo.exe 2>nul

set SRC=dependencies\trading\src\trading_kinds.f90 dependencies\trading\src\trading_strings.f90 dependencies\trading\src\trading_stats.f90 dependencies\trading\src\trading_dependence.f90 dependencies\trading\src\trading_dynamic_beta.f90 dependencies\trading\src\trading_trades.f90 dependencies\trading\src\trading_csa.f90 dependencies\trading\src\trading_curve.f90 dependencies\trading\src\trading_hash_table.f90 dependencies\trading\src\trading_climate.f90 dependencies\trading\src\trading_lottery.f90 dependencies\trading\src\trading_betting.f90 dependencies\trading\src\trading_io.f90 dependencies\trading\src\trading.f90 dependencies\saccr\src\saccr_types.f90 dependencies\saccr\src\saccr_supervisory.f90 dependencies\saccr\src\saccr_core.f90 dependencies\saccr\src\saccr_addon.f90 dependencies\saccr\src\saccr_portfolio.f90 dependencies\saccr\src\saccr_io.f90 dependencies\saccr\src\saccr_examples.f90 dependencies\saccr\src\saccr.f90 src\xva_types.f90 src\xva_math.f90 src\xva_core.f90 src\xva_supervisory.f90 src\xva_exposure.f90 src\xva_regulatory.f90 src\xva_calculator.f90 src\xva.f90

gfortran -std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -J build -I build %SRC% test\test_xva.f90 -o build\test_xva.exe
if errorlevel 1 exit /b 1

gfortran -std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -J build -I build %SRC% app\xva_demo.f90 -o build\xva_demo.exe
if errorlevel 1 exit /b 1

echo Built build\test_xva.exe and build\xva_demo.exe
endlocal
