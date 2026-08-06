# Porting notes

## Laplace objective

For fixed parameters theta, the port minimizes the joint negative log density
`f(h, theta)` over the latent vector `h`. With mode `h_hat` and Hessian `H`, the
integrated negative log likelihood is approximated by

```text
f(h_hat, theta) + 0.5 log(det(H)) - 0.5 n log(2*pi).
```

The AR(1) prior and every observation term involve at most adjacent latent
states, so `H` is tridiagonal. Solves and determinant evaluations are O(n).

## Fixed-parameter transformations

- `sigma_y = exp(log_sigma_y)`
- `sigma_h = exp(log_sigma_h)`
- `phi = tanh(logit_phi/2)`
- `df = exp(log_df_minus_two) + 2`
- `rho = tanh(logit_rho/2)`
- `alpha` is untransformed

## Optimization

A bounded-penalty Nelder-Mead search is used for the 3-5 fixed parameters.
The latent mode is reoptimized at every objective evaluation. This is slower
than TMB automatic differentiation but keeps the port self-contained.

## Uncertainty

The fixed-parameter covariance is the inverse finite-difference Hessian of the
Laplace objective. Natural-scale standard errors use the delta method. Latent
standard errors use the diagonal of the inverse tridiagonal latent Hessian.

## Skew-normal dependency

The upstream R package imports `sn`. The port uses the same centered and
standardized skew-normal parameterization and a hidden-truncation generator.
The relevant formulas were reused under the GPL-3 option of the previous
`sn-fortran` translation.

## Reproducibility

The Fortran RNG is deterministic after `rng%seed(...)`, but it is not R's RNG.
Simulations and optimization results therefore need not match R bit-for-bit.
