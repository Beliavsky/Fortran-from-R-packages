@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build\checked
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set FLAGS=-std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace -fimplicit-none -ffree-line-length-none -O0 -g
set SOURCES=%ROOT%\src\rrcov_kinds.f90 %ROOT%\src\rrcov_types.f90 %ROOT%\src\rrcov_random.f90 %ROOT%\src\rrcov_sort.f90 %ROOT%\src\rrcov_linalg.f90 %ROOT%\src\rrcov_stats.f90 %ROOT%\src\rrcov_robust.f90 %ROOT%\src\rrcov_pca.f90 %ROOT%\src\rrcov_da.f90 %ROOT%\src\rrcov_tests.f90 %ROOT%\src\rrcov_utils.f90 %ROOT%\src\rrcov.f90
gfortran %FLAGS% -c %SOURCES% || exit /b 1
for %%F in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% "%%F" *.o -o "%%~nF.exe" || exit /b 1
  "%%~nF.exe" || exit /b 1
)
endlocal
