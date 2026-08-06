@echo off
setlocal
set FC=gfortran
set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -fcheck=all -fbacktrace
set ROOT=%~dp0
set BUILD=%ROOT%build\windows
set MOD=%BUILD%\mod
set OBJ=%BUILD%\obj
set BIN=%BUILD%\bin

if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%MOD%" "%OBJ%" "%BIN%"

pushd "%OBJ%"
%FC% %FLAGS% -J "%MOD%" -I "%MOD%" -c ^
 "%ROOT%dependencies\qcsis\src\qcsis_kinds.f90" ^
 "%ROOT%dependencies\qcsis\src\qcsis_statistics.f90" ^
 "%ROOT%dependencies\qcsis\src\qcsis.f90" ^
 "%ROOT%dependencies\waveslim\src\waveslim_kinds.f90" ^
 "%ROOT%dependencies\waveslim\src\waveslim_status.f90" ^
 "%ROOT%dependencies\waveslim\src\waveslim_types.f90" ^
 "%ROOT%dependencies\waveslim\src\waveslim_math.f90" ^
 "%ROOT%dependencies\waveslim\src\waveslim_linalg.f90" ^
 "%ROOT%dependencies\waveslim\src\waveslim_filters.f90" ^
 "%ROOT%dependencies\waveslim\src\waveslim_transform_1d.f90" ^
 "%ROOT%dependencies\waveslim\src\waveslim_packet.f90" ^
 "%ROOT%dependencies\waveslim\src\waveslim_transform_nd.f90" ^
 "%ROOT%dependencies\waveslim\src\waveslim_dualtree.f90" ^
 "%ROOT%dependencies\waveslim\src\waveslim_hilbert_stats.f90" ^
 "%ROOT%dependencies\waveslim\src\waveslim_statistics.f90" ^
 "%ROOT%dependencies\waveslim\src\waveslim_denoise.f90" ^
 "%ROOT%dependencies\waveslim\src\waveslim_long_memory.f90" ^
 "%ROOT%dependencies\waveslim\src\waveslim_extended.f90" ^
 "%ROOT%dependencies\waveslim\src\waveslim.f90" ^
 "%ROOT%src\wqc_kinds.f90" "%ROOT%src\wqc_random.f90" ^
 "%ROOT%src\wqc_statistics.f90" "%ROOT%src\wqc.f90"
if errorlevel 1 goto fail

%FC% %FLAGS% -J "%MOD%" -I "%MOD%" "%ROOT%test\test_wqc.f90" *.o -o "%BIN%\test_wqc.exe"
if errorlevel 1 goto fail
"%BIN%\test_wqc.exe"
if errorlevel 1 goto fail

%FC% %FLAGS% -J "%MOD%" -I "%MOD%" "%ROOT%app\wqc_demo.f90" *.o -o "%BIN%\wqc_demo.exe"
if errorlevel 1 goto fail
popd
echo Built %BIN%\wqc_demo.exe
exit /b 0

:fail
popd
echo Build failed.
exit /b 1
