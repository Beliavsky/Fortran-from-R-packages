@echo off
setlocal
call fpm build || exit /b 1
call fpm test || exit /b 1
call fpm run wavethresh_demo || exit /b 1
call fpm run --example basic_dwt || exit /b 1
call fpm run --example denoising || exit /b 1
endlocal
