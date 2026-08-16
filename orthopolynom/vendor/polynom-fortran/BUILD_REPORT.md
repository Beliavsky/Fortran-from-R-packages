# Build report

Compiler used during validation: GNU Fortran 14.2.0.

Validated configurations:

- Checked: Fortran 2018, `-Wall -Wextra -Werror -pedantic`, `-O0`, runtime
  checking, and backtraces.
- Optimized: Fortran 2018, `-Wall -Wextra -Werror -pedantic -O3`.

The package is self-contained and requires no BLAS, LAPACK, C, C++, R, or
external numerical library.
