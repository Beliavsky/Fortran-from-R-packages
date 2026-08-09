# cec2005benchmark-fortran

Modern Fortran translation of the computational code in the R package
`cec2005benchmark` 1.0.4, a wrapper of the CEC 2005 Special Session on
Real-Parameter Optimization benchmark suite.

The package implements all 25 benchmark functions for the four dimensions
supported upstream: 2, 10, 30, and 50.  The official shift vectors, rotation
matrices, composition data, test vectors, and published optima are retained in
`data/`.

## Build

```text
fpm build
fpm test
```

For a strict GNU Fortran build independent of FPM:

```text
scripts/test_gfortran.sh
```

On Windows:

```text
scripts\test_gfortran.bat
```

## Basic use

```fortran
use cec2005benchmark
real(dp) :: x(10), f
integer :: ierr

x = 0.0_dp
f = cec2005_eval(1, x, 'data', .false., ierr)
```

For repeated evaluations, initialize a context once so the benchmark data are
not reread for every point:

```fortran
type(cec2005_context) :: ctx
real(dp) :: x(100,10), f(100)
integer :: ierr

call ctx%init(16, 10, 'data', .false., ierr)
call ctx%evaluate_batch(x, f, ierr)
```

Noise is enabled by default, as in the R package.  Pass `.false.` at
initialization or call `ctx%set_noise(.false.)` for deterministic evaluation.

## Functions

1. Shifted Sphere
2. Shifted Schwefel 1.2
3. Shifted Rotated High-Conditioned Elliptic
4. Shifted Schwefel 1.2 with noise
5. Schwefel 2.6 with optimum on bounds
6. Shifted Rosenbrock
7. Shifted Rotated Griewank
8. Shifted Rotated Ackley with optimum on bounds
9. Shifted Rastrigin
10. Shifted Rotated Rastrigin
11. Shifted Rotated Weierstrass
12. Schwefel 2.13
13. Shifted Expanded Griewank + Rosenbrock (F8F2)
14. Shifted Rotated Expanded Scaffer F6
15. Hybrid Composition Function 1
16. Rotated Hybrid Composition Function 1
17. Rotated Hybrid Composition Function 1 with noise
18. Rotated Hybrid Composition Function 2
19. Rotated Hybrid Composition Function 2 with narrow basin
20. Rotated Hybrid Composition Function 2 with optimum on bounds
21. Rotated Hybrid Composition Function 3
22. Rotated Hybrid Composition Function 3 with high-condition matrix
23. Non-continuous Rotated Hybrid Composition Function 3
24. Rotated Hybrid Composition Function 4
25. Rotated Hybrid Composition Function 4 without bounds

## Validation

The Fortran implementation reproduces all 250 deterministic reference values
in the upstream `test_data_func*.txt` files.  It also reproduces the published
global optimum value of all 25 functions in dimensions 2, 10, 30, and 50.
The worst relative difference observed in those tests with GNU Fortran 14.2.0
was about `2.35e-14`.

See `VALIDATION.md` and `TRANSLATION_COVERAGE.md` for details.

## Licensing

The upstream package declares `GPL (>= 3)`.  This translation is distributed
under GPL-3.0-or-later and retains the complete supplied source tree under
`original/cec2005benchmark-master/`.
