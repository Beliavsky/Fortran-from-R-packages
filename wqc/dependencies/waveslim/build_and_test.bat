@echo off
setlocal
fpm build --profile debug || exit /b 1
fpm test --profile debug || exit /b 1
fpm build --profile release || exit /b 1
fpm test --profile release || exit /b 1
fpm run --profile debug || exit /b 1
fpm run --example image_packet_example --profile debug || exit /b 1
endlocal
