@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build-gfortran-optimized
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\bin"
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -O3

gfortran %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" -c "%ROOT%\src\fastcluster_kinds.f90" -o "%BUILD%\fastcluster_kinds.o" || exit /b 1
gfortran %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" -c "%ROOT%\src\fastcluster_types.f90" -o "%BUILD%\fastcluster_types.o" || exit /b 1
gfortran %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" -c "%ROOT%\src\fastcluster_distances.f90" -o "%BUILD%\fastcluster_distances.o" || exit /b 1
gfortran %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" -c "%ROOT%\src\fastcluster_core.f90" -o "%BUILD%\fastcluster_core.o" || exit /b 1
gfortran %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" -c "%ROOT%\src\fastcluster.f90" -o "%BUILD%\fastcluster.o" || exit /b 1
set OBJECTS="%BUILD%\fastcluster_kinds.o" "%BUILD%\fastcluster_types.o" "%BUILD%\fastcluster_distances.o" "%BUILD%\fastcluster_core.o" "%BUILD%\fastcluster.o"

for %%F in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" "%%F" %OBJECTS% -o "%BUILD%\bin\%%~nF.exe" || exit /b 1
  "%BUILD%\bin\%%~nF.exe" || exit /b 1
)
for %%F in ("%ROOT%\app\*.f90" "%ROOT%\example\*.f90") do (
  gfortran %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" "%%F" %OBJECTS% -o "%BUILD%\bin\%%~nF.exe" || exit /b 1
  "%BUILD%\bin\%%~nF.exe" >nul || exit /b 1
)
endlocal
