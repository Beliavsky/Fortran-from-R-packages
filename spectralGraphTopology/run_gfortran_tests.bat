@echo off
setlocal enabledelayedexpansion
set FC=gfortran
set BUILD=build-gfortran
if exist %BUILD% rmdir /s /q %BUILD%
mkdir %BUILD%\mod %BUILD%\obj %BUILD%\bin
set FLAGS=-std=f2008 -pedantic -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all

for %%F in (src\sgt_kinds.f90 src\sgt_status.f90 src\sgt_types.f90 src\sgt_linalg.f90 src\sgt_operators.f90 src\sgt_utils.f90 src\sgt_updates.f90 src\sgt_objectives.f90 src\sgt_initial_graph.f90 src\sgt_spectral_learning.f90 src\sgt_other_learning.f90 src\spectral_graph_topology.f90) do (
  %FC% %FLAGS% -J %BUILD%\mod -I %BUILD%\mod -c %%F -o %BUILD%\obj\%%~nF.o || exit /b 1
)
set OBJECTS=
for %%F in (%BUILD%\obj\*.o) do set OBJECTS=!OBJECTS! %%F
for %%F in (test\*.f90) do (
  %FC% %FLAGS% -J %BUILD%\mod -I %BUILD%\mod %%F !OBJECTS! -o %BUILD%\bin\%%~nF.exe || exit /b 1
  %BUILD%\bin\%%~nF.exe || exit /b 1
)
%FC% %FLAGS% -J %BUILD%\mod -I %BUILD%\mod app\demo_spectral_graph_topology.f90 !OBJECTS! -o %BUILD%\bin\demo.exe || exit /b 1
%BUILD%\bin\demo.exe || exit /b 1
for %%F in (example\*.f90) do (
  %FC% %FLAGS% -J %BUILD%\mod -I %BUILD%\mod %%F !OBJECTS! -o %BUILD%\bin\%%~nF.exe || exit /b 1
  %BUILD%\bin\%%~nF.exe || exit /b 1
)
endlocal
