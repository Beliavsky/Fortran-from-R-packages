# API coverage

This document records the computational scope of the modern Fortran translation
of R package `wavelets` 0.3-0.2. Machine-readable coverage is in `fpm.toml` and
the summary table is repeated in `README.md`.

## Included computational API

| R function | Fortran mapping | Status | Material differences |
| --- | --- | --- | --- |
| `align` | `wavelets_transform: align` | substantial | Typed transform objects replace R S4 objects and attributes. |
| `dwt` | `wavelets_transform: dwt_named, dwt_with_filter` | substantial | R object/attribute construction and `fast` switch omitted. |
| `dwt.backward` | `wavelets_transform: dwt_backward` | complete | Pure Fortran replaces the upstream C call. |
| `dwt.forward` | `wavelets_transform: dwt_forward` | complete | Pure Fortran replaces the upstream C call. |
| `extend.series` | `wavelets_transform: extend_series` | substantial | Numerical extension modes retained; R class restoration omitted. |
| `idwt` | `wavelets_transform: idwt` | substantial | Numerical inverse retained; R class restoration omitted. |
| `imodwt` | `wavelets_transform: imodwt` | substantial | Numerical inverse retained; R class restoration omitted. |
| `modwt` | `wavelets_transform: modwt_named, modwt_with_filter` | substantial | R object/attribute construction and `fast` switch omitted. |
| `modwt.backward` | `wavelets_transform: modwt_backward` | complete | Pure Fortran replaces the upstream C call. |
| `modwt.forward` | `wavelets_transform: modwt_forward` | complete | Pure Fortran replaces the upstream C call. |
| `mra` | `wavelets_analysis: mra_named, mra_with_filter` | substantial | Numerical MRA retained; R S4/`ts`/data-frame construction omitted. |
| `scalingshift.dwt` | `wavelets_filters: scalingshift_dwt` | complete | Integer phase formulas retained. |
| `waveletshift.dwt` | `wavelets_filters: waveletshift_dwt` | complete | Integer phase formulas retained. |
| `wt.filter` | `wavelets_filters: wt_filter_named, wt_filter_coefficients` | substantial | All named/custom numerical filters retained; R S4 filter class replaced. |
| `wt.filter.equivalent` | `wavelets_filters: wt_filter_equivalent` | complete | Numerical cascade retained. |
| `wt.filter.qmf` | `wavelets_filters: wt_filter_qmf` | complete | QMF convention retained. |
| `wt.filter.shift` | `wavelets_filters: wt_filter_shift` | substantial | Numerical phase rules retained; R polymorphic argument dispatch replaced. |

Coverage status is **substantial**, with **17 of 17 (100%)** exported
computational R functions mapped. The percentage counts mapped computational
entry points and must not be interpreted as 100% R-interface compatibility.

## Excluded from the coverage denominator

The following exported names are deliberately outside computational coverage
because their role is plotting, printing, summaries/formatting, or presentation:

- plotting methods for DWT/MODWT/filter objects;
- printing and summary methods;
- `stackplot`;
- `figure98.wt.filter` and `figure108.wt.filter`;
- `squaredgain.wt.filter`, whose FFT calculation is performed solely to produce
  a filter-response plot and whose externally useful behavior is presentation.

Internal helpers are also excluded from the exported-computational-function
coverage denominator.

## Data model differences

The Fortran API uses `real(dp)` matrices with observations in rows and series in
columns and exposes typed derived types `wt_filter_type`,
`wavelet_transform_type`, and `mra_type`. It does not emulate R S3/S4 dispatch,
`ts` metadata, data-frame column names, arbitrary R attributes, warning objects,
or plotting state. These interface differences are the main reason the package
status remains `substantial` rather than `complete` even though every counted
computational R entry point has a mapping.
