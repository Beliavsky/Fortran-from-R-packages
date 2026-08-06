# Porting notes

## Calibration data

The functions `a`, `a_medium`, and `f` are stored upstream as serialized
`splinefun` closures. The translation extracts each closure's `x`, `y`, `b`,
`c`, and `d` arrays and evaluates the same piecewise cubic polynomial:

```text
y(i) + dx * (b(i) + dx * (c(i) + dx * d(i)))
```

The values are represented by exact IEEE binary64 bit patterns in
`sharpe_rratio_calibration.f90`. `provenance/extract_calibration.py` documents
and reproduces the extraction offsets for the supplied upstream archive.

## Tail fitting

The R routine calls `ghyp::fit.tuv`. The supplied `ghyp` Fortran translation
fits the same skewed Student-t generalized-hyperbolic family by direct maximum
likelihood with Nelder-Mead optimization. The objective is equivalent, but
iteration paths and fitted values need not be bit-for-bit identical to R.

The fitted tail exponent is recovered as `nu = -2*lambda`. If fitting fails,
the estimator falls back to the package's Gaussian sentinel `1.0e13`.

## Random permutations

The original C++ code initializes `std::mt19937` from `std::random_device` and
therefore does not expose reproducibility. The Fortran port uses a local
Park-Miller generator and Fisher-Yates shuffle. An optional `seed` gives fully
reproducible results without modifying Fortran's global random state.

The first permutation remains unshuffled, matching the upstream loop.

## Confidence quantiles

The upstream R function accepts `quantiles=c(0.05,0.95)` but never passes the
values to `computeR0bar`; the C++ defaults 0.025 and 0.975 are therefore always
used. The C++ code also selects order statistics by `floor(p*numPerm)`.

This behavior is retained when `source_compatible` is absent or true. With
`source_compatible=.false.`, `estimate_snr` honors `quantiles` and uses R type-7
interpolation.

## Numerical edge cases

- Empty or wholly non-finite input returns `ok=.false.`.
- Non-finite observations are removed before estimation, as in R.
- A record imbalance with absolute mean below one returns SNR zero.
- Public low-level routines report invalid inputs through result status where
  applicable rather than raising R conditions.

## Combined-work license

The upstream `DESCRIPTION` uses the generic `License: GPL` field. The supplied
`ghyp` dependency is GPL-2.0-or-later. This combined distribution selects
GPL-3.0-only, a version permitted by both components, while retaining the
original license metadata and dependency license texts in the provenance tree.

## One-observation default

For `N=1`, the upstream default `ceiling(3*log(N))` is zero and leads to a
zero-permutation calculation. The Fortran port uses one permutation as the
minimum so the public estimator returns a defined result.
