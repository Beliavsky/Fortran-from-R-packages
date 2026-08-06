# Porting notes

## Array model

R vectors map to `real(dp) :: x(:)`. Multi-asset data map to a matrix with
observations in rows and assets in columns. `estimate_se_matrix` processes each
column independently, matching the R `apply(data, 2, ...)` behavior.

## Tail probabilities

The generic API stores a lower-tail probability in `rpese_options%alpha`.
`var_se` and `es_se` accept confidence levels and perform the same `1-confidence`
conversion as the R wrappers.

## Influence-function dependency

The vendored RPEIF port supplies all influence functions, robust cleaning,
quantiles, partial moments, and robust location estimates. Its
`source_compatibility` option is propagated from RPESE.

## Periodogram fitting

The upstream code excludes the zero Fourier frequency and the Nyquist endpoint,
floors ordinates at `1e-5`, and constructs powers of frequency beginning with an
intercept. This behavior is retained. The zero-frequency fitted spectrum is
obtained by evaluating the fitted log-spectrum polynomial at zero, then divided
by the time-series length to obtain the variance of the estimator.

## Source-compatibility switches

The upstream `DSR` point estimator accepts `rf` but does not subtract it, while
its influence function does. Source-compatible mode preserves that point
estimate; corrected mode uses `(mean-rf)/(sqrt(2)*SemiSD)`.

The upstream adaptive rule assigns full prewhitening weight to negative AR(1)
coefficients because they fall through both nonnegative branches. Source mode
preserves this. Corrected mode clamps the nonnegative coefficient into the
specified transition interval.

## Reproducibility

Cross-validation folds and bootstrap samples use deterministic integer
pseudo-random generators controlled by `rpese_options%seed`.
