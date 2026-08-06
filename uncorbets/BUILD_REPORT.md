# Build report

- Compiler: GNU Fortran 14.2.0
- Language mode: Fortran 2018
- Checked flags: `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`
- Optimized flags: `-O3`
- Warning policy: `-Wall -Wextra -Werror -pedantic`
- External numerical libraries: none
- Runtime dependencies: none

The package is built through the included Makefile. The FPM manifest is also
provided and is checked for unique module/program names and valid TOML syntax.
