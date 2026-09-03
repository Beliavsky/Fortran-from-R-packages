# Build validation

Validation was performed in the supplied execution environment on 2026-09-02.

## Toolchain and native libraries

- GNU Fortran: 14.2.0
- GNU MPFR runtime exercised by the tests: 4.2.2
- GNU GMP runtime: system `libgmp.so.10`
- Fortran standard/options: `-std=f2018 -pedantic`
- Interface checking: `-Wimplicit-interface -Werror=implicit-interface`
- Warnings: `-Wall -Wextra`
- Checked build: `-O0 -fcheck=all`
- Optimized build: `-O2`

The container has the MPFR/GMP runtime libraries but not the unversioned MPFR
development linker name.  For direct validation only, a temporary linker
symlink outside the package tree pointed `libmpfr.so` to the installed
`libmpfr.so.6`.  This validation aid is not included in the package.  On a
normal development installation, FPM uses the native libraries declared by
`link = ["mpfr", "gmp"]` in `fpm.toml`.

## Results

A fresh build directory was used for each configuration.  All package modules,
the deterministic tests, and the example were compiled with explicit module
ordering and linked to the real MPFR/GMP libraries.

- strict checked `-O0` compile: PASS, zero compiler diagnostics
- deterministic tests under `-O0 -fcheck=all`: PASS, zero runtime stderr
- strict optimized `-O2` compile: PASS, zero compiler diagnostics
- deterministic tests under `-O2`: PASS, zero runtime stderr
- `basic_usage` example under `-O2`: PASS, zero compiler diagnostics and zero runtime stderr

The deterministic suite exercises 256-bit construction and rounding, constants,
special functions, exponent controls, exact factorial/binomial results,
probability kernels, summaries, arbitrary-precision matrix products, Romberg
integration, Brent/golden optimization, Brent root finding, qnorm inversion,
and Hooke-Jeeves optimization.  It includes exact/reference checks that exceed
IEEE binary64 precision.

Example output included:

```text
MPFR version: 4.2.2
precision bits: 256
pi: 0.3141592653589793238462643383279502884197169399375105820974944592307816e1
qnorm(0.975): 0.1959963984540054235524594430520551527955550077869548398e1
Bernoulli(10): 0.7575757575757575757575757575757575757576e-1
```

## FPM and fprettify availability

The requested commands were explicitly attempted from the top-level package
directory:

```text
fpm build
fpm test
fpm run --example basic_usage
fpm clean --all
fprettify --version
```

Each returned shell status 127 with `command not found` because this execution
environment contains neither `fpm` nor `fprettify`.  Consequently it would be
incorrect to claim that those literal commands ran successfully here.  The
package was instead validated by the strict direct-gfortran builds described
above, and cleanup was performed manually after the failed `fpm clean --all`
attempt.

## Source and package audits

The maintained Fortran source was checked mechanically for the requested style
and safety properties.  The final audit covers source, tests, and example code.
It verifies:

- free-form source with no line longer than 132 columns;
- no semicolon-separated Fortran statements;
- no `double precision`, `real*8`, `kind(0.0d0)`, or D-exponent literals;
- no self-comparison NaN tests;
- no duplicate Fortran source files;
- every ordinary dummy argument has explicit `INTENT` or `VALUE`;
- dummy arguments are declared one per declaration line and have trailing
  meaningful `!!` FORD comments (procedure-valued dummies are the language-level
  exception to `INTENT`/`VALUE`, but retain their `!!` documentation);
- no copied dependency source or BLAS/LAPACK/ARPACK source;
- no object/module/archive/executable/cache/ZIP products in the package tree.

All ordinary Fortran real variables use the package kind `dp = real64`.  The
three `real(c_double)` declarations in `rmpfr_c_api.f90` are required solely by
the C ABI signatures of `mpfr_set_d`, `mpfr_get_d`, and `mpfr_cmp_d`; they are
not a second numerical working precision.
