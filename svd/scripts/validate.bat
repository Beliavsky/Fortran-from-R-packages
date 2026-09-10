@echo off
setlocal
fpm build || exit /b 1
fpm test || exit /b 1
fpm run svd_demo || exit /b 1
fpm run --example dense_svd || exit /b 1
fpm run --example matrix_free_svd || exit /b 1
python scripts\audit.py || exit /b 1
fpm clean --all || exit /b 1
endlocal
