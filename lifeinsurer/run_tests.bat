@echo off
setlocal
if "%FC%"=="" set FC=gfortran
set ROOT=%~dp0
set BUILD=%ROOT%build-direct
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -O0 -g -J"%BUILD%\mod" -I"%BUILD%\mod"
for %%S in (lifeinsurer_kinds lifeinsurer_types lifeinsurer_helpers lifeinsurer_pv lifeinsurer_cashflows lifeinsurer_actuarial lifeinsurer_profit lifeinsurer_contract lifeinsurer) do %FC% %FLAGS% -c "%ROOT%src\%%S.f90" -o "%BUILD%\obj\%%S.o" || exit /b 1
set OBJS=%BUILD%\obj\lifeinsurer_kinds.o %BUILD%\obj\lifeinsurer_types.o %BUILD%\obj\lifeinsurer_helpers.o %BUILD%\obj\lifeinsurer_pv.o %BUILD%\obj\lifeinsurer_cashflows.o %BUILD%\obj\lifeinsurer_actuarial.o %BUILD%\obj\lifeinsurer_profit.o %BUILD%\obj\lifeinsurer_contract.o %BUILD%\obj\lifeinsurer.o
for %%T in ("%ROOT%test\*.f90") do %FC% %FLAGS% "%%T" %OBJS% -o "%BUILD%\bin\%%~nT.exe" && "%BUILD%\bin\%%~nT.exe" || exit /b 1
endlocal
