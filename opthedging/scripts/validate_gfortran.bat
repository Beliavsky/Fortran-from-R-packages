@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build-validation
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace -ffree-line-length-132 -J"%BUILD%\mod" -I"%BUILD%\mod"

for %%F in (opthedging_kinds opthedging_interpolation opthedging_statistics opthedging_types opthedging_iid opthedging_rng opthedging) do (
  gfortran %FLAGS% -c "%ROOT%\src\%%F.f90" -o "%BUILD%\obj\%%F.o" || exit /b 1
)
set OBJS="%BUILD%\obj\opthedging_kinds.o" "%BUILD%\obj\opthedging_interpolation.o" "%BUILD%\obj\opthedging_statistics.o" "%BUILD%\obj\opthedging_types.o" "%BUILD%\obj\opthedging_iid.o" "%BUILD%\obj\opthedging_rng.o" "%BUILD%\obj\opthedging.o"
for %%F in (test_interpolation test_hedging test_parity test_reference) do (
  gfortran %FLAGS% "%ROOT%\test\%%F.f90" %OBJS% -o "%BUILD%\bin\%%F.exe" || exit /b 1
  "%BUILD%\bin\%%F.exe" || exit /b 1
)
echo validation: PASS
