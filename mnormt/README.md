# mnormt-fortran

**Official CRAN title:** The Multivariate Normal and t Distributions, and Their Truncated Versions

Modern Fortran/FPM translation of the computational part of the R package
`mnormt` 2.1.2.

## Implemented API

- Multivariate normal: `dmnorm`, `dmnorm_many`, `pmnorm`, `rmnorm`, `sadmvn_prob`
- Multivariate t: `dmt`, `dmt_many`, `pmt`, `rmt`, `sadmvt_prob`
- Specialized bivariate/trivariate normal/t probabilities: `biv_nt_prob`, `ptriv_nt`
- Truncated normal: `dmtruncnorm`, `pmtruncnorm`, `rmtruncnorm`
- Truncated t: `dmtrunct`, `pmtrunct`
- Positive-definite inversion: `pd_solve`
- Truncated-normal moment recursion: `recintab`, `mom_mtruncnorm`
- Moment-to-cumulant conversion through order four: `mom2cum`
- Sample Mardia skewness/kurtosis measures: `sample_mardia_measures`

`plot_fxy` and R-specific loading/printing infrastructure are intentionally excluded.

## Build

```sh
fpm test
fpm run --example demo
```

The public API is exposed from module `mnormt` and uses
`dp = kind(1.0d0)`.

## Numerical algorithms

The Genz/Miller/Wichura/Hill probability kernels shipped by upstream `mnormt`
are retained and converted from fixed-form to free-form `.f90`. The modern
module API, linear algebra, RNGs, truncated-normal Gibbs sampler, moment
recursion, and cumulant utilities are independent of R.

The adaptive Genz normal/t rectangle routines retain the upstream maximum
supported dimension of 20. As upstream does for multivariate t rectangle
probabilities, non-integer degrees of freedom are rounded to the nearest
integer in dimensions greater than one.


## FPM compatibility

The retained Genz probability kernels use standard legacy Fortran implicit typing and implicit external procedure interfaces. The manifest explicitly enables these features via `[fortran] implicit-typing = true` and `implicit-external = true`, which is required by current FPM defaults. The modern API modules themselves use explicit typing/interfaces.
