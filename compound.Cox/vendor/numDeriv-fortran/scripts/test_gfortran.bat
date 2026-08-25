@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build\gfortran-strict
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
for %%S in (numderiv_kinds numderiv_types numderiv_callbacks numderiv_core numderiv) do (
  gfortran %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" -c "%ROOT%\src\%%S.f90" -o "%BUILD%\obj\%%S.o" || exit /b 1
)
set OBJECTS="%BUILD%\obj\numderiv_kinds.o" "%BUILD%\obj\numderiv_types.o" "%BUILD%\obj\numderiv_callbacks.o" "%BUILD%\obj\numderiv_core.o" "%BUILD%\obj\numderiv.o"
for %%T in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" "%%T" %OBJECTS% -o "%BUILD%\bin\%%~nT.exe" || exit /b 1
  "%BUILD%\bin\%%~nT.exe" || exit /b 1
)
endlocal
