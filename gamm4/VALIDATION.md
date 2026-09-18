# Validation

Validation date: 2026-09-15.

Compiler available in this environment:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
```

## Dependency situation

`gamm4` is intentionally a sibling-dependent bridge package. Its manifest uses:

```text
mgcv-fortran = { path = "../mgcv" }
lme4 = { path = "../lme4" }
rfortran-linalg = { path = "../rfortran-linalg" }
```

Those sibling checkouts are not mounted in this execution environment. Network
access from the compiler/container environment also cannot resolve GitHub, so
they could not be cloned for a repository-level integration build.

For local compiler validation only, the translation was built against
API-compatible test modules outside the package tree. These test modules expose
the public interfaces used by `gamm4` (`mgcv::smooth_spec_t`, lme4 random-term
and covariance APIs, and the three rfortran-linalg routines). They are **not**
included in the release archive and are not package dependencies.

The actual public dependency signatures were checked against the current
`Fortran-from-R-packages` sources before writing the bridge.

## Strict debug/runtime-check build

The maintained package sources, test, example, and parity driver were compiled
with:

```text
-std=f2018 -O0 -fcheck=all -fbacktrace -Wall -Wextra -Werror -pedantic
```

Results:

```text
All gamm4 tests passed.
```

The example also completed successfully.

## Optimized build

The same source/test/example/parity set was rebuilt with:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -pedantic
```

Results:

```text
All gamm4 tests passed.
```

## Independent NumPy/SciPy parity

`tools/fortran_parity.f90` fits a deterministic Gaussian GAMM containing:

- one parametric intercept;
- one smooth with a one-dimensional null space and two penalized coefficients;
- one ordinary two-level random intercept; and
- ML variance profiling.

`tools/parity_check.py` independently reconstructs and minimizes the same
profiled mixed-model objective with SciPy and separately reconstructs the
upstream `getVb` matrix formula at the Fortran optimum.

Both strict and optimized builds gave essentially the same differences. The
strict-build values were:

```text
objective max/abs difference: 3.476e-08
scale max/abs difference: 5.994e-09
coefficients max/abs difference: 4.183e-07
eta max/abs difference: 9.964e-05
vcov_diag max/abs difference: 1.582e-15
Independent NumPy/SciPy parity checks passed.
```

The optimized `getVb` covariance-diagonal difference was `8.951e-15`.

## Deterministic test coverage

`test/test_gamm4.f90` covers:

- direct `gamm4_get_vb` covariance construction;
- Gaussian GAMM fitting with a smooth and an ordinary grouped random intercept;
- original-basis coefficient/covariance/EDF result assembly;
- Poisson smooth fitting through PIRLS/Laplace machinery; and
- deterministic rejection of unsupported multi-penalty smooths.

## Source-policy audit

`python tools/source_audit.py` reports:

```text
SOURCE AUDIT PASSED: 6 Fortran files; 1/1 coverage mapping consistent.
```

The audit checks maintained Fortran sources for line length, disallowed
semicolon-separated code, legacy real declarations/D exponents, self-comparison
NaN tests, missing trailing dummy-argument FORD comments, duplicate Fortran
content, build products, and manifest/README coverage consistency.

## FPM and fprettify availability

The required FPM commands were explicitly attempted in the package directory:

```text
fpm build
fpm test
fpm clean --all
```

Each returned exit status `127` because `fpm` is not installed on this host.
`fprettify --version` likewise returned `127` because fprettify is not
installed. These commands are therefore **not** reported as successful.

A GitHub connectivity check also failed with:

```text
Could not resolve host: github.com
```

so FPM or the missing sibling packages could not be obtained during this run.
