@echo off
setlocal EnableExtensions EnableDelayedExpansion

rem Always work in the directory containing this script.
cd /d "%~dp0"

set "COMMIT_MESSAGE=Local backup of Fortran packages"

where git >nul 2>&1
if errorlevel 1 (
    echo ERROR: Git was not found in PATH.
    exit /b 1
)

rem Verify that this is the expected parent directory.
for %%P in (
    fGarch-modern-fortran
    longmemo-modern-fortran
    rugarch-modern-fortran
    timsac-modern-fortran
    tseries-modern-fortran
) do (
    if not exist "%%P\fpm.toml" (
        echo ERROR: %%P\fpm.toml was not found.
        echo Run this script from the r_fpm directory.
        exit /b 1
    )
)

rem Add standard generated files and archives to the root ignore file.
call :ensure_ignore "build/"
call :ensure_ignore ".fpm/"
call :ensure_ignore "*.mod"
call :ensure_ignore "*.smod"
call :ensure_ignore "*.o"
call :ensure_ignore "*.obj"
call :ensure_ignore "*.a"
call :ensure_ignore "*.lib"
call :ensure_ignore "*.dll"
call :ensure_ignore "*.so"
call :ensure_ignore "*.dylib"
call :ensure_ignore "*.exe"
call :ensure_ignore "*.out"
call :ensure_ignore "*.zip"

if not exist ".git" (
    echo Initializing one Git repository in:
    echo %CD%
    git init -b main
    if errorlevel 1 exit /b 1
)

git add .
if errorlevel 1 exit /b 1

echo.
echo The following files are staged in the single repository:
echo.
git status --short
echo.
git diff --cached --stat

echo.
set "CONFIRM="
set /p "CONFIRM=Create the local commit shown above? [y/N]: "
if /I not "!CONFIRM!"=="Y" (
    echo No commit was created. Files remain staged for inspection.
    exit /b 0
)

git diff --cached --quiet
if not errorlevel 1 (
    echo There are no new or changed files to commit.
    exit /b 0
)

git commit -m "%COMMIT_MESSAGE%"
if errorlevel 1 (
    echo.
    echo ERROR: The commit failed. If Git requests your identity, run:
    echo   git config --global user.name "Your Name"
    echo   git config --global user.email "you@example.com"
    exit /b 1
)

echo.
echo Local commit completed successfully.
echo No files were uploaded and no network commands were run.
exit /b 0

:ensure_ignore
if not exist ".gitignore" type nul > ".gitignore"
findstr /x /l /c:"%~1" ".gitignore" >nul 2>&1
if errorlevel 1 echo %~1>>".gitignore"
exit /b 0
