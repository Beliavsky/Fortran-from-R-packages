@echo off
rem SPDX-License-Identifier: GPL-2.0-or-later
setlocal
if "%1"=="release" (
  set FLAGS=-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace
) else (
  set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace
)
if not exist build mkdir build
set SRC=src\fbonds_kinds.f90 src\fbonds_linalg.f90 src\fbonds_csv.f90 src\fbonds_term_structure.f90
gfortran %FLAGS% -Jbuild -Ibuild -c %SRC%
if errorlevel 1 exit /b 1
gfortran %FLAGS% -Jbuild -Ibuild test\test_term_structure.f90 fbonds_kinds.o fbonds_linalg.o fbonds_csv.o fbonds_term_structure.o -llapack -lblas -o build\test_term_structure.exe
if errorlevel 1 exit /b 1
build\test_term_structure.exe
if errorlevel 1 exit /b 1
gfortran %FLAGS% -Jbuild -Ibuild app\demo_fbonds.f90 fbonds_kinds.o fbonds_linalg.o fbonds_csv.o fbonds_term_structure.o -llapack -lblas -o build\demo_fbonds.exe
gfortran %FLAGS% -Jbuild -Ibuild app\fit_csv.f90 fbonds_kinds.o fbonds_linalg.o fbonds_csv.o fbonds_term_structure.o -llapack -lblas -o build\fit_csv.exe
gfortran %FLAGS% -Jbuild -Ibuild example\fit_example.f90 fbonds_kinds.o fbonds_linalg.o fbonds_csv.o fbonds_term_structure.o -llapack -lblas -o build\fit_example.exe
if errorlevel 1 exit /b 1
build\demo_fbonds.exe > nul
build\fit_csv.exe data\example_yield.csv ns > nul
build\fit_csv.exe data\example_yield.csv svensson l1 > nul
build\fit_example.exe > nul
echo build and tests passed.
