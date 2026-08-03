@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set SOURCES=jrvfinance_kinds.f90 jrvfinance_types.f90 jrvfinance_dates.f90 jrvfinance_roots.f90 jrvfinance_cashflows.f90 jrvfinance_bonds.f90 jrvfinance_options.f90 jrvfinance.f90
call :run debug "-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow" || exit /b 1
call :run release "-std=f2018 -O3 -Wall -Wextra -Wpedantic -Werror" || exit /b 1
echo All jrvFinance-fortran validations passed.
exit /b 0
:run
set NAME=%~1
set FLAGS=%~2
set BUILD=%ROOT%\build-%NAME%
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod"
pushd "%BUILD%"
for %%S in (%SOURCES%) do gfortran %FLAGS% -J mod -I mod -c "%ROOT%\src\%%S" || (popd & exit /b 1)
for %%T in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% -J mod -I mod "%%T" *.o -o "%%~nT.exe" || (popd & exit /b 1)
  "%%~nT.exe" || (popd & exit /b 1)
)
for %%T in ("%ROOT%\example\*.f90" "%ROOT%\app\*.f90") do (
  gfortran %FLAGS% -J mod -I mod "%%T" *.o -o "%%~nT.exe" || (popd & exit /b 1)
  "%%~nT.exe" >nul || (popd & exit /b 1)
)
popd
exit /b 0
