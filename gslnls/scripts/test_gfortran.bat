@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build-strict
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod"
mkdir "%BUILD%\bin"
cd /d "%BUILD%"
if "%FC%"=="" set FC=gfortran
if "%FFLAGS%"=="" set FFLAGS=-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all

%FC% %FFLAGS% -Jmod -Imod -c "%ROOT%\src\gslnls_kinds.f90" -o gslnls_kinds.o || exit /b 1
%FC% %FFLAGS% -Jmod -Imod -c "%ROOT%\src\gslnls_types.f90" -o gslnls_types.o || exit /b 1
%FC% %FFLAGS% -Jmod -Imod -c "%ROOT%\src\gslnls_linalg.f90" -o gslnls_linalg.o || exit /b 1
%FC% %FFLAGS% -Jmod -Imod -c "%ROOT%\src\gslnls_loss.f90" -o gslnls_loss.o || exit /b 1
%FC% %FFLAGS% -Jmod -Imod -c "%ROOT%\src\gslnls_core.f90" -o gslnls_core.o || exit /b 1
%FC% %FFLAGS% -Jmod -Imod -c "%ROOT%\src\gslnls_multistart.f90" -o gslnls_multistart.o || exit /b 1
%FC% %FFLAGS% -Jmod -Imod -c "%ROOT%\src\gslnls_stats.f90" -o gslnls_stats.o || exit /b 1
%FC% %FFLAGS% -Jmod -Imod -c "%ROOT%\src\gslnls.f90" -o gslnls.o || exit /b 1
set OBJS=gslnls_kinds.o gslnls_types.o gslnls_linalg.o gslnls_loss.o gslnls_core.o gslnls_multistart.o gslnls_stats.o gslnls.o
for %%F in ("%ROOT%\test\*.f90") do (
  echo -- %%~nxF
  %FC% %FFLAGS% -Jmod -Imod "%%F" %OBJS% -o "bin\%%~nF.exe" || exit /b 1
  "bin\%%~nF.exe" || exit /b 1
)
for %%F in ("%ROOT%\example\*.f90") do (
  echo -- %%~nxF
  %FC% %FFLAGS% -Jmod -Imod "%%F" %OBJS% -o "bin\%%~nF.exe" || exit /b 1
  "bin\%%~nF.exe" || exit /b 1
)
echo All strict tests and examples passed.
endlocal
