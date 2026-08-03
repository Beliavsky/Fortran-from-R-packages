# API

All public names are available through `use gnorm`.

## Kinds and status codes

- `dp`: double-precision real kind
- `i8`: 64-bit integer kind
- `gnorm_success`
- `gnorm_invalid_argument`
- `gnorm_numerical_failure`

## Distribution functions

```fortran
value = dgnorm(x, mu, alpha, beta, log_density)
value = pgnorm(q, mu, alpha, beta, lower_tail, log_probability)
value = qgnorm(p, mu, alpha, beta, lower_tail, log_probability)
```

All arguments after the first are optional. Defaults are `mu=0`, `alpha=1`,
`beta=1`, lower tail, and ordinary rather than logarithmic output.

The functions are elemental and accept scalar or conformable array arguments.
Invalid parameters or probabilities return IEEE NaN. Quantiles at probabilities
zero and one return negative and positive infinity.

## Random generation

```fortran
values = rgnorm(n, mu, alpha, beta, seed, status)
call rgnorm_fill(values, mu, alpha, beta, seed, status)
```

`seed` is an optional 64-bit integer. Supplying it makes the generated sequence
reproducible. Without it, the generator is initialized from `system_clock`.

## Moments

```fortran
mean_value = gnorm_mean(mu)
variance = gnorm_variance(alpha, beta)
```

The variance is
`alpha**2 * Gamma(3/beta) / Gamma(1/beta)`.

## Numerical helpers

The following are public for testing and reuse:

- `regularized_gamma_p`
- `regularized_gamma_q`
- `inverse_regularized_gamma_p`
- `gnorm_rng_state`
