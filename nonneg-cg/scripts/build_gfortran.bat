@echo off
setlocal
if not exist build\manual mkdir build\manual
gfortran -std=f2018 -O2 -Wall -Wextra -Wpedantic -Wimplicit-interface -Werror=implicit-interface -J build\manual -I build\manual -c src\nonneg_cg.f90 -o build\manual\nonneg_cg.o || exit /b 1
gfortran -std=f2018 -O2 -Wall -Wextra -Wpedantic -Wimplicit-interface -Werror=implicit-interface -J build\manual -I build\manual build\manual\nonneg_cg.o test\test_nonneg_cg.f90 -o build\manual\test_nonneg_cg.exe || exit /b 1
gfortran -std=f2018 -O2 -Wall -Wextra -Wpedantic -Wimplicit-interface -Werror=implicit-interface -J build\manual -I build\manual build\manual\nonneg_cg.o example\rosenbrock.f90 -o build\manual\rosenbrock.exe || exit /b 1
build\manual\test_nonneg_cg.exe || exit /b 1
build\manual\rosenbrock.exe || exit /b 1
