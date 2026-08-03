@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0
set BUILD=%ROOT%build-gfortran
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
if "%FC%"=="" set FC=gfortran
if "%FFLAGS%"=="" set FFLAGS=-std=f2018 -Wall -Wextra -Werror -pedantic -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -O0 -g
set SOURCES=fints_kinds fints_status fints_types fints_linalg fints_special fints_summary fints_time_series fints_finance fints_dates fints_apca fints_arma fints_arima fints
set OBJECTS=
for %%S in (%SOURCES%) do (
  %FC% %FFLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" -c "%ROOT%src\%%S.f90" -o "%BUILD%\obj\%%S.o" || exit /b 1
  set OBJECTS=!OBJECTS! "%BUILD%\obj\%%S.o"
)
for %%T in ("%ROOT%test\*.f90") do (
  %FC% %FFLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" "%%T" !OBJECTS! -o "%BUILD%\bin\%%~nT.exe" || exit /b 1
  "%BUILD%\bin\%%~nT.exe" || exit /b 1
)
for %%T in ("%ROOT%app\*.f90" "%ROOT%example\*.f90") do (
  %FC% %FFLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" "%%T" !OBJECTS! -o "%BUILD%\bin\%%~nT.exe" || exit /b 1
  "%BUILD%\bin\%%~nT.exe" >nul || exit /b 1
)
echo All FinTS-fortran tests and runnable targets passed.
