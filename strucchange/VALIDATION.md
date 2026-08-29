# Validation

Validation was performed on 2026-08-29 with GNU Fortran 14.2.0.

## Source and API checks

The release tree was checked for:

- free-form `.f90` maintained Fortran sources only;
- a single real kind, `dp` from `r_kinds`, in `src/`, `test/`, and `example/`;
- no `double precision`, `real*8`, `kind(0.0d0)`, or `d0`/`D` exponent literals;
- no semicolon-separated statements;
- no source line longer than 132 characters;
- no duplicate Fortran source content;
- no copied `rfortran-core`, `rfortran-linalg`, BLAS, LAPACK, ARPACK,
  `r.f90`, or `r_mod.f90` source;
- no object/module/executable/cache/build/ZIP artifacts inside the package.

The package's uses of `r_kinds`, `r_distributions`, and `r_linalg` were checked
against the current shared-module APIs in the target repository. In particular,
`dp` is `real64`, and `rfortran-linalg` uses the same `real64` kind for its
LAPACK-backed interfaces.

## Compiler validation

Because the execution environment did not provide FPM, an API-compatible
external harness for the two sibling dependencies was used only outside the
release directory to perform a strict independent compilation of every package
source file. The harness is not included in the ZIP.

Compiler flags used for package/test/example compilation were:

```text
-std=f2018 -Wall -Wextra -Wimplicit-interface -fcheck=all
```

The validation linked against the system LAPACK and BLAS libraries. Every
`strucchange` source compiled without warnings. The only compiler warnings were
unused-local warnings in the external temporary linear-algebra harness.

The deterministic test executable reported:

```text
All strucchange tests passed.
```

The examples reported:

```text
BIC-selected number of breaks: 1
Breakpoints: 20
Maximum absolute OLS-CUSUM:   1.3514E+00
```

Tests include independently calculated reference values for recursive
residuals, supF/aveF/expF statistics, one/two/three-break RSS values, segmented
regression coefficients, a breakpoint confidence interval, OLS-CUSUM and
OLS-MOSUM values, generalized score-process standardization, supLM/maxMOSUM and
categorical-L2 functionals, asymptotic p-values, and monitoring critical values.

## FPM limitation in this environment

The requested commands are the intended release gate:

```sh
fpm build
fpm test
fpm clean --all
```

However, no `fpm` executable is installed in this sandbox. A download attempt
for the current standalone FPM binary also failed because executable downloads
are unavailable from the runtime. Therefore these three commands could not be
truthfully recorded as successful here.

`tools/release_check.py` performs those exact FPM commands together with the
source/dependency/artifact checks and should be run after placing `strucchange/`
beside `rfortran-core/` and `rfortran-linalg/` in the target repository.
