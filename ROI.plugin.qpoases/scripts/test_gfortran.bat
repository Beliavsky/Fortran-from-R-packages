@echo off
setlocal
if "%FC%"=="" set FC=gfortran
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
set ROOT=%~dp0..
set BUILD=%ROOT%\build-strict
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"

%FC% %FLAGS% -J"%BUILD%\mod" -c "%ROOT%\src\qpoases_kinds.f90" -o "%BUILD%\obj\qpoases_kinds.o" || exit /b 1
%FC% %FLAGS% -I"%BUILD%\mod" -J"%BUILD%\mod" -c "%ROOT%\src\qpoases_types.f90" -o "%BUILD%\obj\qpoases_types.o" || exit /b 1
%FC% %FLAGS% -I"%BUILD%\mod" -J"%BUILD%\mod" -c "%ROOT%\src\qpoases_linalg.f90" -o "%BUILD%\obj\qpoases_linalg.o" || exit /b 1
%FC% %FLAGS% -I"%BUILD%\mod" -J"%BUILD%\mod" -c "%ROOT%\src\qpoases_active_set.f90" -o "%BUILD%\obj\qpoases_active_set.o" || exit /b 1
%FC% %FLAGS% -I"%BUILD%\mod" -J"%BUILD%\mod" -c "%ROOT%\src\qpoases_solver.f90" -o "%BUILD%\obj\qpoases_solver.o" || exit /b 1
%FC% %FLAGS% -I"%BUILD%\mod" -J"%BUILD%\mod" -c "%ROOT%\src\qpoases.f90" -o "%BUILD%\obj\qpoases.o" || exit /b 1
%FC% %FLAGS% -I"%BUILD%\mod" -J"%BUILD%\mod" -c "%ROOT%\src\roi_qpoases.f90" -o "%BUILD%\obj\roi_qpoases.o" || exit /b 1

set OBJS="%BUILD%\obj\qpoases_kinds.o" "%BUILD%\obj\qpoases_types.o" "%BUILD%\obj\qpoases_linalg.o" "%BUILD%\obj\qpoases_active_set.o" "%BUILD%\obj\qpoases_solver.o" "%BUILD%\obj\qpoases.o" "%BUILD%\obj\roi_qpoases.o"

for %%F in ("%ROOT%\test\*.f90") do (
  %FC% %FLAGS% -I"%BUILD%\mod" "%%F" %OBJS% -o "%BUILD%\bin\%%~nF.exe" || exit /b 1
  "%BUILD%\bin\%%~nF.exe" || exit /b 1
)

for %%F in ("%ROOT%\example\*.f90") do (
  %FC% %FLAGS% -I"%BUILD%\mod" "%%F" %OBJS% -o "%BUILD%\bin\%%~nF.exe" || exit /b 1
  "%BUILD%\bin\%%~nF.exe" || exit /b 1
)

echo All strict tests and examples passed.
