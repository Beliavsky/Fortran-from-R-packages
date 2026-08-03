@echo off
setlocal enabledelayedexpansion
set FC=gfortran
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -Wno-compare-reals -Wno-intrinsic-shadow -fcheck=all -fbacktrace
set BUILD=build\checked
if exist %BUILD% rmdir /s /q %BUILD%
mkdir %BUILD%

for %%M in (metrics_kinds metrics_utils metrics_regression metrics_binary_classification metrics_classification metrics_information_retrieval metrics_time_series metrics metrics_test_support) do (
  %FC% %FLAGS% -J %BUILD% -I %BUILD% -c src\%%M.f90 -o %BUILD%\%%M.o || exit /b 1
)
set OBJS=%BUILD%\metrics_kinds.o %BUILD%\metrics_utils.o %BUILD%\metrics_regression.o %BUILD%\metrics_binary_classification.o %BUILD%\metrics_classification.o %BUILD%\metrics_information_retrieval.o %BUILD%\metrics_time_series.o %BUILD%\metrics.o %BUILD%\metrics_test_support.o
for %%F in (test\*.f90) do (
  %FC% %FLAGS% -J %BUILD% -I %BUILD% %%F %OBJS% -o %BUILD%\%%~nF.exe || exit /b 1
  %BUILD%\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90 app\*.f90) do (
  %FC% %FLAGS% -J %BUILD% -I %BUILD% %%F %OBJS% -o %BUILD%\%%~nF.exe || exit /b 1
  %BUILD%\%%~nF.exe || exit /b 1
)
