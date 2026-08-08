@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build-strict
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -J"%BUILD%\mod" -I"%BUILD%\mod"
for %%S in (rceim_kinds rceim_random rceim_utils rceim_benchmarks rceim) do (
  gfortran %FLAGS% -c "%ROOT%\src\%%S.f90" -o "%BUILD%\obj\%%S.o" || exit /b 1
)
set OBJS="%BUILD%\obj\rceim_kinds.o" "%BUILD%\obj\rceim_random.o" "%BUILD%\obj\rceim_utils.o" "%BUILD%\obj\rceim_benchmarks.o" "%BUILD%\obj\rceim.o"
for %%T in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% "%%T" %OBJS% -o "%BUILD%\bin\%%~nT.exe" || exit /b 1
  "%BUILD%\bin\%%~nT.exe" || exit /b 1
)
for %%E in ("%ROOT%\example\*.f90") do (
  gfortran %FLAGS% "%%E" %OBJS% -o "%BUILD%\bin\%%~nE.exe" || exit /b 1
  "%BUILD%\bin\%%~nE.exe" || exit /b 1
)
echo All strict tests and examples passed.
endlocal
