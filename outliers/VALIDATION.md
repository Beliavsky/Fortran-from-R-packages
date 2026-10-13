# Validation

## Reference environment

Reference calculations used actual R 4.5.0 and the installed CRAN outliers 0.15
package from the pinned official archive. The deterministic regeneration program
is `tools/generate_references.R`. The installed R package is required only to
regenerate fixtures; ordinary FPM builds and tests do not require R or Python.

The checked-in fixture set contains:

* 22,126 distribution records: all usable Dixon n/type combinations, automatic
  type selection, all probability-table nodes, interpolation, both tails of
  extrapolation, and reverse queries; all Grubbs variants up to n=1000 where
  supported; Cochran n=2, 2.5, 3, 5, 10, 30, 100 and k=2, 3, 5, 10, 30.
* 576 individual-test configurations, including all type/opposite/two-sided
  combinations, supplied/default chi-square variance, both endpoint directions,
  pairs, tied values, NaNs, singleton/two-point samples and constant samples.
* 1,080 score configurations across five score types, raw/CDF/threshold modes,
  IQR limits, missing values, sample-size boundaries and degenerate scales.

The FPM test program performs **32,433 checks** with zero failures. Additional
checks cover independent polynomial interpolation, exact table nodes, forward/
reverse inverses where mathematically valid, tied endpoint regression, table
anomalies, invalid sizes, removal/replacement, matrix shapes, masks, ragged columns,
Cochran sample counts and grouped-data interfaces.

## Tolerances and explicit differences

Ordinary analytic calculations and scores use absolute tolerances 2e-11 or
smaller. Dixon forward interpolation uses 2e-10; reverse interpolation uses 5e-8
because direct Lagrange evaluation and R lm/poly can differ near ill-conditioned
anomalous table stencils. Independent review observed a maximum finite reverse
Dixon difference of 2.38e-8, with forward differences about 1e-13.

Grubbs type 11 is inverted algebraically. R uses uniroot with its default root
accuracy; a tolerance of 8e-5 for the CDF, or 1.6e-4 for doubled test p-values,
reflects that source approximation rather than a looser tolerance for the shared
t distribution. Forward critical values are checked at 2e-11. Independent review
observed maximum valid-root CDF disagreement about 3.04e-5.

The source type-11 clipping/bracket discontinuity has an unavoidable branch
ambiguity when the statistic differs from the threshold by one backend-rounding
ulp. The implementation recognizes equality within 128 machine epsilons times
statistic scale; its exact lower-boundary result is NaN. At this precise boundary,
tests explicitly accept the reference's zero or bracket-error NaN category and
verify the implementation returns zero/NaN. This does not change any tabular
sample-size limit or the original probability thresholds. Interior cases and
ordinary out-of-domain cases are tested strictly.

Frozen R fixtures are left unchanged for invalid Grubbs samples. Tests explicitly
check the documented Fortran policy instead: n<3 and constant samples must return
status 1 and NaNs, rather than preserving R's finite G without usable probability
or its accidental constant-data p=0. NaN/infinity results in score fixtures are
compared by category, not discarded.

## Independent review

A separate reviewer compared 4,508 Dixon values, thousands of Grubbs/Cochran
values, a 558-row Grubbs inverse sweep, 529 API cases and six Cochran matrix/group
examples against actual R. It identified and verified fixes for local duplicate
Dixon coordinates and two-point t scores. No remaining numerical/API blockers
were reported. Cochran matrix, unequal-group, inlying/outlying and tied-extrema
comparisons agreed within 1e-14.

## Build matrix and source audit

Executed with GNU Fortran and FPM on Linux:

```sh
fpm build
fpm test
fpm run --example demo
fpm build --profile debug --flag "-O0 -g -Wall -Wextra -fcheck=all -fbacktrace"
fpm test --profile debug --flag "-O0 -g -Wall -Wextra -fcheck=all -fbacktrace"
fpm run --profile debug --flag "-O0 -g -Wall -Wextra -fcheck=all -fbacktrace" --example demo
fpm build --profile release
fpm test --profile release
fpm run --profile release --example demo
```

All nine stages pass. Windows portability is supported by standard free-form
Fortran, relative sibling FPM dependencies and no system numerical libraries;
Windows itself was not executed in this validation environment.

The audit verifies preserved upstream files against all 27 MD5 manifest entries;
coverage metadata agrees with 14 exact exported R names; each maintained dummy
has its own INTENT declaration and meaningful trailing FORD comment; all real
arithmetic uses the imported dp; no semicolon-separated statements or
self-comparison NaN tests occur; there are no duplicate Fortran sources, vendored
dependencies, or BLAS/LAPACK links. Existing sibling dependency code is not
rewritten or bundled. Build products are removed with `fpm clean --all` after
verification.
