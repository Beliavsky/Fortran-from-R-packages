# Design notes

## FFT

R's `fft()` calls were replaced by a native complex FFT. Power-of-two transforms use an iterative radix-2 algorithm. Other lengths use Bluestein convolution, whose internal convolution length is padded to a power of two. Forward transforms are unnormalized; `fft_inverse_raw` matches R's unnormalized inverse transform convention.

## Gaussian simulation

`sim_gauss` preserves the original circulant-embedding construction. Small negative eigenvalues attributable to floating-point roundoff are clipped to zero; materially negative values cause an error.

## Spectra

The formulas in `spec_fgn`, `b_spec_fgn`, and `spec_arima` follow the R implementation directly. The normalized fGn spectrum is divided by the geometric-mean scale `theta1` exactly as in the package.

## Covariance matrices

`ceta_fgn` and `ceta_arima` use the package's forward finite differences of the log spectrum. The derivative cross-product is inverted through a Cholesky factorization.

## Whittle estimation

The one-parameter fGn problem uses a bounded Brent-style scalar minimizer over `[0.1,0.99]`. The multi-parameter fARIMA problem uses unconstrained Nelder-Mead, matching the absence of parameter constraints in the R call to `optim()`.

## FEXP estimation

For a Gamma GLM with log link, Fisher-scoring weights are constant. The implementation performs IRLS through repeated least-squares updates and estimates dispersion from Pearson residuals. Student-t p-values use a regularized incomplete-beta implementation.

The polynomial columns are produced with modified Gram-Schmidt orthogonalization. This spans the same model space as R's `poly(ffr,j)` and reproduces the fitted spectrum and H estimate closely. Nuisance-polynomial coefficient scaling is basis-dependent.
