@echo off
setlocal
if exist build rmdir /s /q build
if exist objects rmdir /s /q objects
mkdir build
mkdir objects
set FFLAGS=-std=f2008 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -O0
gfortran %FFLAGS% -J build -I build -c src\mnb_kinds.f90 src\mnb_types.f90 src\mnb_math.f90 src\mnb_optimizer.f90 src\mnb_core.f90 src\mnb_simulation.f90 src\mnb_residuals.f90 src\mnb_influence.f90 src\mnb_envelope.f90 src\mnb.f90
move /y *.o objects\ >nul
for %%F in (test\*.f90) do (
  gfortran %FFLAGS% -J build -I build objects\*.o %%F -o build\%%~nF.exe || exit /b 1
  build\%%~nF.exe || exit /b 1
)
gfortran %FFLAGS% -J build -I build objects\*.o example\demo_mnb.f90 -o build\demo_mnb.exe || exit /b 1
build\demo_mnb.exe || exit /b 1
endlocal
