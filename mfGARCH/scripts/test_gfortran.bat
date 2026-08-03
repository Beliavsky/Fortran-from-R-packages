@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build\gfortran-debug
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\obj" || exit /b 1
mkdir "%BUILD%\bin" || exit /b 1
cd /d "%BUILD%\obj" || exit /b 1
set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -fbacktrace

gfortran %FLAGS% -c ^
 "%ROOT%\src\mfgarch_kinds.f90" ^
 "%ROOT%\src\mfgarch_status.f90" ^
 "%ROOT%\src\mfgarch_types.f90" ^
 "%ROOT%\src\mfgarch_math.f90" ^
 "%ROOT%\src\mfgarch_components.f90" ^
 "%ROOT%\src\mfgarch_low_level.f90" ^
 "%ROOT%\src\mfgarch_optimization.f90" ^
 "%ROOT%\src\mfgarch_random.f90" ^
 "%ROOT%\src\mfgarch_simulation.f90" ^
 "%ROOT%\src\mfgarch_prediction.f90" ^
 "%ROOT%\src\mfgarch_fit.f90" ^
 "%ROOT%\src\mfgarch_io.f90" ^
 "%ROOT%\src\mfgarch.f90" || exit /b 1

for %%F in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% *.o "%%F" -o "%BUILD%\bin\%%~nF.exe" || exit /b 1
  "%BUILD%\bin\%%~nF.exe" || exit /b 1
)
for %%F in ("%ROOT%\example\*.f90" "%ROOT%\app\*.f90") do (
  gfortran %FLAGS% *.o "%%F" -o "%BUILD%\bin\%%~nF.exe" || exit /b 1
  pushd "%BUILD%"
  "%BUILD%\bin\%%~nF.exe" || exit /b 1
  popd
)
endlocal
