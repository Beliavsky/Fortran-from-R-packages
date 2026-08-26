@echo off
rem Configure FPM/GFortran to use Windows BLAS, LAPACK, and ARPACK import libraries.
rem Run this once in a cmd.exe session before building packages that link them.

if defined R_FPM_LINALG_ROOT if exist "%R_FPM_LINALG_ROOT%\lib\liblapack.dll.a" goto configure

set "R_FPM_LINALG_ROOT="
for /f "delims=" %%E in ('where octave-cli.exe 2^>nul') do call :try_octave_exe "%%E"
if defined R_FPM_LINALG_ROOT goto configure

for /d %%D in ("C:\Programs\Octave-*\mingw64") do if exist "%%~fD\lib\liblapack.dll.a" set "R_FPM_LINALG_ROOT=%%~fD"
if defined R_FPM_LINALG_ROOT goto configure

for /d %%D in ("C:\Program Files\GNU Octave\Octave-*\mingw64") do if exist "%%~fD\lib\liblapack.dll.a" set "R_FPM_LINALG_ROOT=%%~fD"
if defined R_FPM_LINALG_ROOT goto configure

for /d %%D in ("C:\Programs\JAGS\JAGS-*\x64") do if exist "%%~fD\lib\liblapack.dll.a" set "R_FPM_LINALG_ROOT=%%~fD"
if defined R_FPM_LINALG_ROOT goto configure

echo ERROR: No compatible BLAS/LAPACK import libraries were found.
echo Install GNU Octave for Windows, or set R_FPM_LINALG_ROOT to a directory containing bin and lib.
exit /b 1

:try_octave_exe
for %%D in ("%~dp1..") do if exist "%%~fD\lib\liblapack.dll.a" set "R_FPM_LINALG_ROOT=%%~fD"
exit /b 0

:configure
set "R_FPM_LINALG_FLAG=-L%R_FPM_LINALG_ROOT:\=/%/lib"
if not defined FPM_LDFLAGS goto set_link_flag
echo(%FPM_LDFLAGS%| %SystemRoot%\System32\findstr.exe /l /c:"%R_FPM_LINALG_FLAG%" >nul
if errorlevel 1 set "FPM_LDFLAGS=%R_FPM_LINALG_FLAG% %FPM_LDFLAGS%"
goto set_path

:set_link_flag
set "FPM_LDFLAGS=%R_FPM_LINALG_FLAG%"

:set_path
echo(;%PATH%;| %SystemRoot%\System32\findstr.exe /i /l /c:";%R_FPM_LINALG_ROOT%\bin;" >nul
if errorlevel 1 set "PATH=%R_FPM_LINALG_ROOT%\bin;%PATH%"
echo FPM will use numerical libraries from "%R_FPM_LINALG_ROOT%".
echo FPM_LDFLAGS=%FPM_LDFLAGS%
exit /b 0
