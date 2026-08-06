# QCSIS modern Fortran

A self-contained modern Fortran translation of the computational code in the
R package **QCSIS 0.1**. It implements quantile correlation, composite
quantile correlation, QC-SIS, and CQC-SIS. The project uses the Fortran Package
Manager (FPM), has no external dependencies, and contains no plotting code.

## Implemented API

The public module is `qcsis_mod`.

- `qc(x, y, tau [, stat, errmsg])` returns `type(qc_result)` with `tau` and
  `rho` arrays.
- `cqc(x, y [, stat, errmsg])` returns the scalar composite quantile
  correlation.
- `qcsis(x, y, d [, stat, errmsg])` uses the R default quantile grid
  `1:(n-1)/n`.
- `qcsis(x, y, tau, d [, stat, errmsg])` uses a caller-supplied quantile grid.
- `cqcsis(x, y, d [, stat, errmsg])` performs composite quantile-correlation
  sure independence screening.

Screening functions return `type(screening_result)` with:

- `w`: one screening weight per predictor;
- `selected`: the one-based indices of the top `d` predictors.

The translation reproduces R's default type-7 sample quantile, R's sample
standard deviation convention (denominator `n - 1`), and the original strict
comparison `y - quantile < 0`.

## Build and test

```text
fpm build
fpm test
fpm run --example basic_qcsis
```

A direct GNU Fortran build is also possible:

```text
gfortran -std=f2018 -Wall -Wextra -Wpedantic -fcheck=all \
  src/qcsis_kinds.f90 src/qcsis_statistics.f90 src/qcsis.f90 \
  test/test_qcsis.f90 -o test_qcsis
./test_qcsis
```

## Input requirements

- `x` and `y` must be finite.
- Vector inputs must have equal length and at least two observations.
- Every predictor must have a positive finite sample standard deviation.
- Each probability in `tau` must be strictly between zero and one.
- `d` must be between one and the number of predictor columns.

When supplied, `stat` is zero on success. On failure it is nonzero and
`errmsg` describes the invalid input. Failed derived-type results contain
zero-length arrays.

## Upstream and license

The original R package was written by Xuejun Ma, Jingxiao Zhang, and Jingke
Zhou and released under GPL-2. This translation preserves that license and
attribution. The upstream `DESCRIPTION`, `NAMESPACE`, and R computational
sources are retained under `original/` for provenance.

See `LICENSE` for the GNU General Public License version 2.
