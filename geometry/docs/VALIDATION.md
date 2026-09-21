# Validation

The maintained Fortran sources are designed for Fortran 2018 and GNU Fortran on Windows or Unix-like systems.

Validation in the translation environment used GNU Fortran 14.2.0. The release checks compile with both a strict runtime-checking configuration and an optimized configuration, then run the two deterministic test programs and the example.

Strict configuration used:

```text
gfortran -std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace
```

Optimized configuration used:

```text
gfortran -std=f2018 -O2
```

The literal commands requested for FPM are also attempted before packaging:

```text
fpm build
fpm test
fpm run --example geometry_example
fpm clean --all
```

The current execution sandbox does not provide an `fpm` executable, so those literal commands return exit status 127 (`fpm: command not found`). This is an environment limitation, not reported as an FPM pass. The same source files, tests, and example are compiled and run directly with GNU Fortran instead.

The computational-geometry implementations deliberately avoid the upstream Qhull C code. The replacement hull and Delaunay algorithms are combinatorial and are intended for small and moderate problem sizes. They are not a performance-compatible replacement for Qhull on thousands of points or high dimensions. See `docs/API_COVERAGE.md` for function-level compatibility notes.
