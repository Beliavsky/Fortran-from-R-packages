# wavelets

Modern free-form Fortran translation of the computational core of the R package
`wavelets` 0.3-0.2 by Eric Aldrich. The upstream package computes discrete
wavelet transforms (DWT), maximal-overlap discrete wavelet transforms (MODWT),
their inverses, filters, phase shifts, series extensions, alignment, and
multiresolution analyses.

This translation is intended to live as the top-level `wavelets` directory in
<https://github.com/Beliavsky/Fortran-from-R-packages>.

## Build and test

A Fortran Package Manager (FPM) installation and a Fortran 2018-capable compiler
such as gfortran are sufficient. No separately installed BLAS or LAPACK library
is required.

```text
fpm build
fpm test
fpm run --example basic_wavelets
```

The maintained Fortran source is free form. The project has no external FPM
dependencies and does not vendor BLAS, LAPACK, ARPACK, R compatibility modules,
or translated R packages.

## Public API

The umbrella module is `wavelets` and re-exports the common real kind `dp`, the
public result/filter types, and the translated numerical procedures.

- Filters: `wt_filter`, `wt_filter_named`, `wt_filter_coefficients`,
  `wt_filter_qmf`, `wt_filter_equivalent`, `wt_filter_shift`,
  `waveletshift_dwt`, `scalingshift_dwt`.
- DWT: `dwt`, `dwt_named`, `dwt_with_filter`, `dwt_forward`, `dwt_backward`,
  `idwt`.
- MODWT: `modwt`, `modwt_named`, `modwt_with_filter`, `modwt_forward`,
  `modwt_backward`, `imodwt`.
- Utilities and analysis: `extend_series`, `align`, `mra`, `mra_named`,
  `mra_with_filter`.

High-level series inputs are `real(dp)` matrices with observations in rows and
series in columns. Typed Fortran objects replace the upstream R S4 objects and
R-specific attributes. Status-returning APIs use optional/integer `ierr`
arguments rather than R warnings/errors.

All 25 named filters implemented by the upstream package are included: Haar
(`haar`/`d2`), Daubechies `d4` through `d20`, least-asymmetric `la8` through
`la20`, best-localized `bl14`, `bl18`, `bl20`, and Coiflets `c6`, `c12`, `c18`,
`c24`, `c30`.

## Numerical compatibility notes

The four upstream C transform kernels are translated into pure Fortran rather
than retained as C interoperability wrappers. Periodic and reflection boundary
handling, default DWT level selection, MODWT normalization for named filters,
filter phase conventions, alignment, and multiresolution reconstruction follow
the upstream algorithms.

The R `fast` argument is not exposed because the Fortran implementation always
uses compiled numerical kernels. R S3/S4 dispatch, `ts`/data-frame attributes,
plotting, printing, summaries, and presentation-only filter figures are outside
the translation scope. Inverse high-level transforms retain the upstream
five-decimal rounding behavior. See `API_COVERAGE.md` for details.

## Translation coverage

Package status: **substantial**.

**17 of 17 (100%)** exported computational R functions have a meaningful
Fortran mapping. This fraction measures mapped computational functions, not
complete R interface or object-system compatibility. Plotting, printing,
formatting, interactive/display, presentation-only functions, and internal
helpers are excluded from the denominator under the project's stated coverage
convention.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `align` | `wavelets_transform` | `align` | substantial |
| `dwt` | `wavelets_transform` | `dwt_named`, `dwt_with_filter` | substantial |
| `dwt.backward` | `wavelets_transform` | `dwt_backward` | complete |
| `dwt.forward` | `wavelets_transform` | `dwt_forward` | complete |
| `extend.series` | `wavelets_transform` | `extend_series` | substantial |
| `idwt` | `wavelets_transform` | `idwt` | substantial |
| `imodwt` | `wavelets_transform` | `imodwt` | substantial |
| `modwt` | `wavelets_transform` | `modwt_named`, `modwt_with_filter` | substantial |
| `modwt.backward` | `wavelets_transform` | `modwt_backward` | complete |
| `modwt.forward` | `wavelets_transform` | `modwt_forward` | complete |
| `mra` | `wavelets_analysis` | `mra_named`, `mra_with_filter` | substantial |
| `scalingshift.dwt` | `wavelets_filters` | `scalingshift_dwt` | complete |
| `waveletshift.dwt` | `wavelets_filters` | `waveletshift_dwt` | complete |
| `wt.filter` | `wavelets_filters` | `wt_filter_named`, `wt_filter_coefficients` | substantial |
| `wt.filter.equivalent` | `wavelets_filters` | `wt_filter_equivalent` | complete |
| `wt.filter.qmf` | `wavelets_filters` | `wt_filter_qmf` | complete |
| `wt.filter.shift` | `wavelets_filters` | `wt_filter_shift` | substantial |

The machine-readable mapping is in `[extra.translation]` and
`[[extra.translation.function]]` entries in `fpm.toml`; those entries and this
table are intended to remain synchronized.

## Tests

The deterministic test program checks:

- filter construction, QMF conventions, equivalent filters, and phase shifts;
- all series-extension modes;
- the low-level Haar DWT and MODWT kernels against fixed reference values;
- DWT/IDWT round trips with periodic and reflection boundaries;
- MODWT/IMODWT round trips with periodic and reflection boundaries;
- alignment and inverse alignment; and
- DWT- and MODWT-based multiresolution reconstruction identities.

## License and provenance

The upstream package declares `GPL (>= 2)`. This translation is distributed as
GPL-2.0-or-later to match that license. See `COPYING`, `LICENSE`, `NOTICE`, and
the retained computational source snapshots under `upstream/`.

The original upstream author and maintainer is Eric Aldrich
(`<ealdrich@gmail.com>`). The upstream CRAN package metadata records version
0.3-0.2 and publication on 2020-02-17. No claim is made that this Fortran port is
an official release of the upstream author.
