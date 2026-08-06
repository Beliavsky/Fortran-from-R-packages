# wqc modern Fortran

A modern Fortran/FPM translation of the computational code in the R package
`wqc` 0.1.2, **Wavelet Quantile Correlation Analysis**.

The translation implements:

- `quantile_correlation_analysis` for one pair of time series;
- `apply_quantile_correlation` for a matrix whose first column is the
  reference series and whose remaining columns are targets;
- `quantile_correlation`, a vector-quantile wrapper around the translated
  `QCSIS::qc` routine;
- MODWT multiresolution analysis through the attached modern Fortran
  translation of `waveslim`;
- Gaussian-surrogate 95% confidence intervals;
- R-compatible type-7 quantiles and sample standard deviations;
- optional deterministic random seeding, validation, tests, and an example.

Plotting code is not translated.

## Build

The dependency translations are bundled as local FPM path dependencies, so no
network access is required.

```text
fpm test
fpm run
```

With GNU Fortran but without FPM:

```text
./build_and_test.sh
```

On Windows with GNU Fortran available on `PATH`:

```text
build_and_test.bat
```

## Pairwise API

```fortran
use wqc, only : dp, quantile_correlation_analysis, wqc_pair_result

type(wqc_pair_result) :: fit
integer :: stat
character(len=:), allocatable :: errmsg

fit = quantile_correlation_analysis(x, y, [0.1_dp, 0.5_dp, 0.9_dp], &
   wf='la8', j_levels=4, n_sim=1000, seed=12345, &
   stat=stat, errmsg=errmsg)
```

The result arrays have shape `(levels, number_of_quantiles)`:

- `fit%estimated_qc`
- `fit%ci_lower`
- `fit%ci_upper`

The quantile vector is stored in `fit%quantiles`.

## Multi-series API

```fortran
use wqc, only : apply_quantile_correlation, dp, wqc_multi_result

type(wqc_multi_result) :: fit
real(dp) :: data(n, p)

fit = apply_quantile_correlation(data, [0.1_dp, 0.5_dp, 0.9_dp], &
   j_levels=4, n_sim=500, seed=12345)
```

`data(:,1)` is the reference series. Each `fit%series(i)` corresponds to
`data(:,i+1)`. Optional `series_names` may contain one name per target column.

## Computational equivalence

The wavelet decomposition, quantile-correlation formula, sample standard
deviation, and confidence-interval quantiles follow the R implementation.
The R package uses R's global normal generator; this port uses a local
Park-Miller generator with Box-Muller normal deviates. Therefore seeded Monte
Carlo confidence intervals are reproducible within this port but are not
expected to match R bit for bit. Estimated quantile correlations do not depend
on the random generator.

## Licenses

The top-level wqc translation preserves the upstream GPL-3 license. The
bundled dependencies remain separate FPM packages and retain their original
licenses:

- `dependencies/qcsis`: GPL-2
- `dependencies/waveslim`: BSD-3-Clause

See `LICENSE`, `NOTICE`, and the license files in each dependency directory.

License compatibility note: upstream wqc declares GPL-3 while importing the
GPL-2 QCSIS package. This port mirrors that dependency relationship and
keeps the source packages separate. Redistribution of a linked binary
containing both components may require a separate license review.
