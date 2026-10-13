# Validation

## Reproducible setup

Verified on 2026-10-11 with GNU Fortran 14.2.0 and FPM 0.13.0 (Linux).
The upstream reference package was built and loaded from CRAN mirror commit
`2321fc1788c387ad06e531d408c31b5d840b2286`, splines2 0.5.4, using R 4.5.0,
Rcpp 1.0.14 and RcppArmadillo 14.4.2-1. These R dependencies are used only
for reference generation and are not included in this Fortran package.
The shared repository dependencies were inspected at commit
`cbd437d6d92f5481706fd8ea5a7829aaff7b6e24`.

Run from this package directory:

```sh
fpm build
fpm test
fpm run --example
fpm build --profile debug
fpm test --profile debug
fpm run --example --profile debug
fpm build --profile release
fpm test --profile release
fpm run --example --profile release
```

All builds, tests and examples passed. Additional debug builds/tests/examples
with `--flag '-fcheck=all -fbacktrace -Wall -Wextra'` also passed. No unsafe
floating-point optimization flags are used. Windows execution was not available;
the maintained source and manifests are platform-neutral free-form Fortran and
sibling FPM paths, without external BLAS/LAPACK links.

## Numerical checks

`test/test_reference.f90` reads the committed deterministic numeric
`test/reference.dat`. It requires no R installation. There are **1,190 upstream
reference cases**, containing complete numerical matrices:

- B/M/I/C degrees 0 through 5, arbitrary derivative orders including those beyond
  polynomial degree, complete/reduced bases, integration and C scaling.
- Interior points, knots, exact boundaries and polynomial extrapolation.
- Periodic B/M degrees 0 through 5, derivatives and integrated bases across
  negative and positive periods, unequal intervals and non-unit boundaries.
- Natural and cardinal cubic bases with zero through five interior knots,
  derivatives zero through five, integrals and linear/quadratic extrapolation.
- Bernstein degrees zero through ten, derivatives, integrals and exact endpoints.
- Repeated B/M/I/C knots through multiplicity degree+1.

Acceptance: maximum elementwise `abs(actual-expected)/max(1,abs(expected))`
below `2e-10`. Observed worst error across default, debug and release profiles:
**3.6379788070917130e-12**. For the known upstream reduced scaled C-derivative
failure, the generator uses upstream's full derivative then removes column one;
this deliberate upstream-bug repair is recorded in `API_COVERAGE.md`.

`test/test_identities.f90` runs **231 independent identity/API assertions**:
partition of unity and zero derivative sums, B nonnegativity, unit M integrals,
I monotonicity, C convexity and normalization, derivative/integral finite
differences, I/M/C calculus relations, periodicity and accumulated areas,
natural endpoint constraints and linear extrapolation, nsk cardinality,
Bernstein/B equivalence, aliases and wrapper defaults, retained prediction
settings, NaN handling, type-7/trim/fallback knot selection, repeated-boundary
sorting, and invalid-argument statuses. It also verifies the mathematically
valid no-interior-knot periodic degree-zero/one extension for which upstream's
df constructor may fail.

An independent review additionally probed **768 automatic-knot/df/trim cases**
against the same compiled R package. All numerically valid upstream results
matched within `2e-10`. It identified the periodic-M explicit `df >= degree`
guard, which was then added, and the intentionally retained no-knot periodic
upstream-error repair documented above.

## Regenerating the fixtures

Install the pinned upstream splines2 0.5.4 (with Rcpp/RcppArmadillo), then run:

```sh
Rscript tools/generate_reference.R
fpm test
```

The generator checks the installed package version; the precise pinned source
commit is recorded in its fixture header and `NOTICE.md`. Reference files are
plain text, column-major, with family/settings and shape information preceding
each matrix. Production builds and tests never invoke R.

## Source and packaging checks

The maintained Fortran consists of two modules, two tests and one example.
The static audit found no duplicate Fortran source, forbidden semicolon-separated
statements, self-comparison NaN tests, missing dummy INTENT/VALUE or meaningful
FORD comments, inconsistent real kinds, or coverage-count inconsistencies.
`dp` comes only from `r_kinds`; all dependencies are sibling paths.
No dependency source, Rcpp, Armadillo, BLAS or LAPACK is copied here.

Before release packaging, after independent review and any reruns, execute
`fpm clean --all` and ensure no build products, cache directories, ZIP files,
executables, objects or module files remain in the package directory.
