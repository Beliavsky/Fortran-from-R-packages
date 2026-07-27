@echo off
setlocal
set FC=gfortran
set FLAGS=-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace
if not exist build\windows mkdir build\windows
%FC% %FLAGS% -Jbuild\windows -c src\backtest_kinds.f90 -o build\windows\backtest_kinds.o || exit /b 1
%FC% %FLAGS% -Jbuild\windows -Ibuild\windows -c src\backtest_math.f90 -o build\windows\backtest_math.o || exit /b 1
%FC% %FLAGS% -Jbuild\windows -Ibuild\windows -c src\backtest_bucket.f90 -o build\windows\backtest_bucket.o || exit /b 1
%FC% %FLAGS% -Jbuild\windows -Ibuild\windows -c src\backtest_portfolio.f90 -o build\windows\backtest_portfolio.o || exit /b 1
%FC% %FLAGS% -Jbuild\windows -Ibuild\windows -c src\backtest_engine.f90 -o build\windows\backtest_engine.o || exit /b 1
%FC% %FLAGS% -Jbuild\windows -Ibuild\windows -c src\backtest_summary.f90 -o build\windows\backtest_summary.o || exit /b 1
set OBJS=build\windows\backtest_kinds.o build\windows\backtest_math.o build\windows\backtest_bucket.o build\windows\backtest_portfolio.o build\windows\backtest_engine.o build\windows\backtest_summary.o
%FC% %FLAGS% -Ibuild\windows app\demo_backtest.f90 %OBJS% -o build\windows\demo_backtest.exe || exit /b 1
%FC% %FLAGS% -Ibuild\windows app\backtest_csv.f90 %OBJS% -o build\windows\backtest_csv.exe || exit /b 1
%FC% %FLAGS% -Ibuild\windows example\overlap_example.f90 %OBJS% -o build\windows\overlap_example.exe || exit /b 1
echo Windows release build passed.
