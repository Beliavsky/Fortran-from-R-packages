# Validation record

Date: 2026-10-11. Platform: Linux x86-64. Compiler: GNU Fortran 14.2.0.
Build tool: FPM 0.13.0. Reference runtime: R 4.5.0, with posterior 1.7.0
installed from the pinned source in NOTICE.md. Native Windows execution was
not available; the build uses standard free-form Fortran, sibling FPM paths
and no system numerical-library links.

## Final independent build matrix

All nine build/test/example stages passed after the final source correction:

| Profile | Build | Tests | Examples |
|---|---|---|---|
| Default FPM | Passed | Passed | Passed |
| Checked debug | Passed | Passed | Passed |
| Release | Passed | Passed | Passed |

Commands from the package directory:

```sh
fpm build
fpm test
fpm run --example
fpm build --profile debug --flag '-g -O0 -fcheck=all -fbacktrace -Wall -Wextra'
fpm test --profile debug --flag '-g -O0 -fcheck=all -fbacktrace -Wall -Wextra'
fpm run --example --profile debug --flag '-g -O0 -fcheck=all -fbacktrace -Wall -Wextra'
fpm build --profile release
fpm test --profile release
fpm run --example --profile release
```

The independent matrix cleaned the build tree before each profile. An
additional implementation-time checked run passed with `-std=f2018`,
`-Wimplicit-interface` and `-Werror=line-truncation`. No fast-math or
finite-math-only option was used.

## Tests and R references

Reported numerical checks per full Fortran test run:

- Diagnostics: **411** comparisons/properties, including actual pinned-R
  fixtures for independent, positively autocorrelated, antithetic, shifted,
  scaled, tied, odd/short, one-chain, constant, missing and infinite draws.
- Draw/rvar computations: **159** deterministic assertions, including sorted
  duplicate subsets, ties-to-even thinning, supplied-uniform resampling,
  broadcasting, log weights, empty reductions, missing/infinite quantiles and
  numeric/logical reductions.
- Pareto/predictive computations: **2745** numerical comparisons/properties,
  plus status/invalid-input assertions. The extra invalid long-tail-name
  regression is an uncounted guard assertion, not part of that printed total.
  Cases cover bounded/exponential/heavy/tied and unsorted observations,
  zero-weight fitting, four smoothing modes, weighted PIT, all three uniformity
  tests, inverse-distribution identities and discrete summaries.
- Numeric transformations: callback ordering/mutation, matrix-margin
  reductions, diagonal construction/recycling, matching and integer/logical
  coercion assertions. These are reported without a numerical count.

Thus **3315 counted** checks pass, in addition to the transformation and guard
assertions. This is a test count, not a count of unique R functions.

The four supplied R scripts execute successfully against pinned posterior:

```sh
Rscript test/reference_diagnostics.R
Rscript test/reference_pareto.R
Rscript test/reference_draws.R
Rscript test/reference_transforms.R
```

The first two regenerate checked-in numeric fixtures. The draw and transform
scripts separately validate their matching deterministic cases with **60**
and **12** actual-R comparisons. R is not required for `fpm test`.

A separate review compared **2533 diagnostic values** across 149 additional
cases, and **640 Pareto/PIT values** across 40 weighted/unweighted cases. Those
checks passed at stated tolerances without changing the release fixtures.

## Resolved optimized-build discrepancy

The first release pass exposed scalar/vector EXP rounding at a tied log-tail
cutoff. Scalar and vectorized EXP differed by one ULP, producing a tiny
negative excess even though sorted tail observations mathematically imply a
nonnegative excess. The validated log-tail path now enforces that nonnegative
domain before GPD fitting. The original R smoothing references then pass in
all profiles at the original tolerances. The issue was independently
reproduced before correction; it was not hidden by tolerance changes.

## Focused floating-point guard test

A focused `is_constant` probe passed with
`-g -O0 -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`.
It verifies repeated infinity, repeated NaN, a finite range larger than the
largest representable real, and NaN/infinite/negative/zero tolerances. Explicit
separate finiteness and sign guards avoid assuming short-circuit evaluation.
The corresponding seven regression assertions are included in the diagnostic
test program. This is a focused guard check, not a claim that every deliberately
undefined numerical summary can run with all IEEE exceptions trapped.

## Static and coverage checks

The independent static audit passed with no errors or warnings:

- **11** maintained free-form Fortran files
- **287** dummy arguments checked for explicit INTENT and meaningful trailing
  FORD documentation; procedure-pointer callbacks also declare INTENT
- One named real kind, imported from `r_kinds`
- No disallowed semicolon-separated statements or self-comparison NaN tests
- No duplicate maintained Fortran sources or vendored dependency source
- No system BLAS/LAPACK links or embedded compatibility modules
- Coverage reconciles to **91/113**, with exactly **22** untranslated names
  and individually documented **27** exclusions from 140 explicit exports
- README, API coverage and TOML mappings agree; mapped procedures are public
- Original R provenance sources, metadata, license and SHA-256 manifest retained

Release cleanup uses `fpm clean --all` after final review. Build products,
executables, object/module files, caches and dependency sources are excluded
from the deliverable.
