@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build-gfortran-debug
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod"
mkdir "%BUILD%\obj"
set FC=gfortran
set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -J"%BUILD%\mod" -I"%BUILD%\mod"
call :compile "%ROOT%\dependencies\lbfgs\src\lbfgs_kinds.f90" lbfgs_kinds || exit /b 1
call :compile "%ROOT%\dependencies\lbfgs\src\lbfgs_status.f90" lbfgs_status || exit /b 1
call :compile "%ROOT%\dependencies\lbfgs\src\lbfgs_solver.f90" lbfgs_solver || exit /b 1
call :compile "%ROOT%\dependencies\lbfgs\src\lbfgs.f90" lbfgs || exit /b 1
call :compile "%ROOT%\dependencies\lbfgsb3\src\lbfgsb3_core.f90" lbfgsb3_core || exit /b 1
call :compile "%ROOT%\dependencies\lbfgsb3\src\lbfgsb3.f90" lbfgsb3 || exit /b 1
call :compile "%ROOT%\src\rcppnumerical_kinds.f90" rcppnumerical_kinds || exit /b 1
call :compile "%ROOT%\src\rcppnumerical_callbacks.f90" rcppnumerical_callbacks || exit /b 1
call :compile "%ROOT%\src\rcppnumerical_gk_data.f90" rcppnumerical_gk_data || exit /b 1
call :compile "%ROOT%\src\rcppnumerical_integration_1d.f90" rcppnumerical_integration_1d || exit /b 1
call :compile "%ROOT%\src\rcppnumerical_cuhre.f90" rcppnumerical_cuhre || exit /b 1
call :compile "%ROOT%\src\rcppnumerical_optimization.f90" rcppnumerical_optimization || exit /b 1
call :compile "%ROOT%\src\rcppnumerical_logistic.f90" rcppnumerical_logistic || exit /b 1
call :compile "%ROOT%\src\rcppnumerical.f90" rcppnumerical || exit /b 1
set OBJS="%BUILD%\obj\lbfgs_kinds.o" "%BUILD%\obj\lbfgs_status.o" "%BUILD%\obj\lbfgs_solver.o" "%BUILD%\obj\lbfgs.o" "%BUILD%\obj\lbfgsb3_core.o" "%BUILD%\obj\lbfgsb3.o" "%BUILD%\obj\rcppnumerical_kinds.o" "%BUILD%\obj\rcppnumerical_callbacks.o" "%BUILD%\obj\rcppnumerical_gk_data.o" "%BUILD%\obj\rcppnumerical_integration_1d.o" "%BUILD%\obj\rcppnumerical_cuhre.o" "%BUILD%\obj\rcppnumerical_optimization.o" "%BUILD%\obj\rcppnumerical_logistic.o" "%BUILD%\obj\rcppnumerical.o"
%FC% %FLAGS% "%ROOT%\test\test_rcppnumerical.f90" %OBJS% -o "%BUILD%\test_rcppnumerical.exe" || exit /b 1
"%BUILD%\test_rcppnumerical.exe" || exit /b 1
%FC% %FLAGS% "%ROOT%\example\integration_optimization.f90" %OBJS% -o "%BUILD%\integration_optimization.exe" || exit /b 1
"%BUILD%\integration_optimization.exe" || exit /b 1
%FC% %FLAGS% "%ROOT%\example\fast_lr_example.f90" %OBJS% -o "%BUILD%\fast_lr_example.exe" || exit /b 1
"%BUILD%\fast_lr_example.exe" || exit /b 1
exit /b 0

:compile
%FC% %FLAGS% -c %1 -o "%BUILD%\obj\%2.o"
exit /b %errorlevel%
