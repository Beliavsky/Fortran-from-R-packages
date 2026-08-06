# Porting notes

## Architecture

The R package combines R S3 classes, C++ kernels, univariate `tsgarch` models,
optimization packages, and presentation methods.  The Fortran port separates
this into numerical modules:

- `tsmarch_types`: specifications and result containers
- `tsmarch_linalg`: self-contained matrix and random-number support
- `tsmarch_dcc`: DCC filtering, likelihood, estimation, simulation, forecasting
- `tsmarch_copula`: Gaussian and Student copula-GARCH
- `tsmarch_ica`: whitening, FastICA, and RADICAL-style ICA
- `tsmarch_gogarch`: GO-GARCH estimation, moments, and forecasting
- `tsmarch_utils`: covariance, distribution, diagnostic, and matrix utilities
- `tsmarch_risk`: portfolio simulation risk measures
- `tsmarch_fft`: portable convolution and distribution evaluation
- `tsmarch`: umbrella module

The required `ghyp`, `tsdistributions`, and `tsgarch` computational modules are
vendored in `src` so the package builds without external Fortran dependencies.

## DCC parameterization

DCC coefficients are represented directly by `dcc_parameters`.  Estimation
maps them to an unconstrained vector and back, enforcing nonnegative dynamic
coefficients and a stable aggregate persistence.  Asymmetric DCC uses the
negative-shock cross product and its sample target matrix.

The unconditional correlation target can be based on Pearson, Spearman, or a
Kendall-derived correlation, selected by `dcc_spec%correlation_method`.

## Marginal models

Each column is fitted independently through the bundled `tsgarch` numerical
API.  The DCC and copula layers consume standardized residuals and conditional
standard deviations from those fits.  This replaces the R package's
`tsgarch.multi_estimate` and S3 conversion machinery.

## Copula implementation

Marginal residuals are transformed through their fitted CDFs.  Gaussian copulas
use normal quantiles; Student copulas use standardized Student quantiles and a
shape parameter.  The correlation recursion is shared with the DCC engine.

## ICA and GO-GARCH

Whitening uses a symmetric eigendecomposition.  FastICA supports symmetric and
deflation updates.  The RADICAL entry point performs deterministic pairwise
rotations minimizing a spacing-entropy criterion.  Since ICA is identifiable
only up to sign, scale, and permutation, component ordering need not match R.

Factor GARCH models are fitted independently with the bundled `tsgarch` API,
and the mixing matrix maps factor variances and higher moments back to asset
space.

## Convolution

`fft_convolution` exposes the same computational role as the upstream GO-GARCH
FFT distribution methods.  For maximum portability and deterministic tests it
uses a zero-padded direct discrete Fourier transform.  This is accurate for the
small and medium grids targeted by the tests, but applications with very large
grids may prefer replacing the internal DFT with FFTW or another optimized FFT.

## Numerical inference

DCC inference is computed with finite-difference Hessians and observation score
vectors.  The result type stores the Hessian, its inverse covariance estimate,
standard errors, and score matrix.  These are numerical equivalents of R's
`numDeriv` and sandwich-method paths, not automatic differentiation.

## Random numbers

The Fortran intrinsic generator is seeded deterministically through an expanded
64-bit seed.  Results are reproducible within this implementation but do not
share R's Mersenne-Twister stream.
