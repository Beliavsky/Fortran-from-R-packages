@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build-strict
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
set FLAGS=-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
for %%M in (lsei_kinds lsei_types lsei_linalg lsei_nnls lsei_solver lsei_utils lsei) do (
  gfortran %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" -c "%ROOT%\src\%%M.f90" -o "%BUILD%\obj\%%M.o" || exit /b 1
)
for %%T in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% -I "%BUILD%\mod" "%%T" "%BUILD%\obj\lsei_kinds.o" "%BUILD%\obj\lsei_types.o" "%BUILD%\obj\lsei_linalg.o" "%BUILD%\obj\lsei_nnls.o" "%BUILD%\obj\lsei_solver.o" "%BUILD%\obj\lsei_utils.o" "%BUILD%\obj\lsei.o" -o "%BUILD%\bin\%%~nT.exe" || exit /b 1
  "%BUILD%\bin\%%~nT.exe" || exit /b 1
)
for %%T in ("%ROOT%\example\*.f90") do (
  gfortran %FLAGS% -I "%BUILD%\mod" "%%T" "%BUILD%\obj\lsei_kinds.o" "%BUILD%\obj\lsei_types.o" "%BUILD%\obj\lsei_linalg.o" "%BUILD%\obj\lsei_nnls.o" "%BUILD%\obj\lsei_solver.o" "%BUILD%\obj\lsei_utils.o" "%BUILD%\obj\lsei.o" -o "%BUILD%\bin\%%~nT.exe" || exit /b 1
  "%BUILD%\bin\%%~nT.exe" || exit /b 1
)
endlocal
