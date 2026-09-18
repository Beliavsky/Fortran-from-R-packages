# Validation

Validation was performed in the translation environment with:

- GNU Fortran (Debian 14.2.0-19) 14.2.0
- Python 3 with NumPy and SciPy for independent numerical references

## Shared `mclust` dependency

`InitClust` now declares the sibling FPM dependency:

```text
mclust-fortran = { path = "../mclust" }
```

The repository's current `mclust_hierarchical` public interface was audited before integration.  It exports
`hc_result`, `hc_fit`, and `hclass`, and its `dp` kind is binary64-compatible with `otrimle`'s `real64` kind.
The actual sibling checkout is not mounted in this execution environment and FPM is unavailable, so direct gfortran
validation used a non-shipped interface test double with the same procedure signatures.  A second non-shipped test
double forces `hc_fit` failure and verifies the Manhattan average-link fallback.  Neither test double is included in
the release archive.

## Direct strict builds

The maintained `otrimle` sources plus the non-shipped `mclust` interface test double were compiled in dependency
order without system BLAS/LAPACK libraries.

Debug/runtime-check configuration:

```text
gfortran -std=f2018 -O0 -fcheck=all -fbacktrace -Wall -Wextra -Werror -pedantic
```

Result:

```text
All otrimle tests passed.
```

The example also ran successfully. Its representative output was:

```text
status code: 2
improper log likelihood:  -6.891644
component proportions:   0.100023   0.499988   0.399989
cluster 1 mean:    -2.0000    -2.0000
cluster 2 mean:     2.0250     2.0500
```

Optimized configuration:

```text
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -pedantic
```

The complete deterministic test suite, example, and independent parity suite also passed in this configuration.

A separate debug/runtime-check build with `hc_fit` forced to fail passed the complete deterministic tests, exercising
the upstream Manhattan average-link fallback path.

## Compatibility regression tests

The deterministic suite now contains fixed numerical checks for:

- the upstream default `robustbase::scaleTau2(..., mu.too=TRUE)` location and scale;
- the R `stats::density()` weighted-binning/Gaussian-convolution/interpolation path used by `kerndensmeasure`; and
- both successful hierarchical-initialization dispatch and the Manhattan fallback integration.

For the fixed `scaleTau2` sample, an independent implementation gives the same robust location and scale to roughly
machine precision.  For the fixed density sample, all 100 Fortran density ordinates agree with an independent
reconstruction of R's binned-convolution algorithm to a maximum absolute difference of about `2.2e-16`.

## Independent numerical parity

`tools/parity_driver.f90` exports a fitted model, density diagnostic, and `scaleTau2` reference values.
`tools/check_parity.py` recomputes the corresponding quantities independently with NumPy/SciPy.
Both strict builds passed with the following observed differences:

```text
improper log likelihood:       7.461e-14
posterior probability:         4.258e-11
weighted chi-square criterion: 2.498e-16
R-compatible density measure:  0.000e+00
scaleTau2 robust location:      5.551e-17
scaleTau2 robust scale:         1.110e-15
```

The independent check also verifies the requested global covariance eigenvalue-ratio bound and maximum mean
noise-posterior constraint.

## Static policy audit

Run:

```text
python tools/source_audit.py
```

Current result:

```text
SOURCE AUDIT PASSED: 11 Fortran files; 12/12 coverage mapping consistent.
```

The audit checks, among other things:

- the single shared `dp` definition in maintained `otrimle` source;
- free-form line length at or below 132 columns;
- no executable semicolon-separated statements;
- no `double precision`, `real*8`, `kind(0.0d0)`, or D-exponent literals;
- no self-comparison NaN idioms;
- explicit `INTENT`/`VALUE`, one dummy per declaration, and trailing `!!` FORD comments;
- no duplicate Fortran files;
- no copied `r.f90` or `r_mod.f90`;
- no compiler/build products or nested ZIP files in the package tree; and
- consistency of translation counts in `fpm.toml` and `README.md`.

## FPM and fprettify availability

The required FPM commands were explicitly attempted from the package root:

```text
fpm build
fpm test
fpm clean --all
```

All three returned shell exit status 127 because `fpm` is not installed in this execution environment.  Therefore
this document does **not** claim an FPM run that did not occur.  The direct gfortran builds above validate the
maintained package sources and the exact `mclust` call signatures; a target-repository FPM build additionally
requires the sibling `../mclust` package, as intended by the manifest.

`fprettify --version` likewise returned exit status 127 because `fprettify` is not installed. Source formatting was
checked by the static policy audit instead.

## Release archive check

The release package tree contains one top-level directory named `otrimle`, no copied `mclust`, `robustbase`, BLAS,
or LAPACK source, and no compiler products.  The final ZIP is re-extracted into a fresh directory and the maintained
sources are recompiled with the non-shipped `mclust` interface test double in both strict configurations before the
archive is released.
