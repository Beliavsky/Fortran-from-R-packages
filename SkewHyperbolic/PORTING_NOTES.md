# Porting notes

- Plotting, Q-Q/P-P/tail-plot helpers, R S3 methods, datasets, histogram bookkeeping, and formula/presentation code are intentionally omitted.
- The density follows the upstream Bessel-K formula. The `beta = 0` branch is the corresponding scaled Student-t distribution.
- `pskewhyp` integrates the density after a finite-interval transformation; `qskewhyp` uses bracket expansion plus bisection. The upstream package offers spline acceleration as an alternative; that is a performance convenience rather than a distinct numerical model and is not reproduced.
- `rskewhyp` directly ports the defining mixture: `W = 1/G`, with `G ~ Gamma(nu/2, scale=2/delta^2)`, and `X = mu + beta W + sqrt(W) Z`.
- Integer moments are evaluated analytically from that mixture, rather than through `DistributionUtils::momRecursion` and `gammaRawMom`.
- `skewhypFit` is translated as `skewhyp_fit`. It uses the same transformed parameter space `(mu, log(delta), beta, log(nu))`, but replaces R's `optim`/`nlm` choices with a deterministic derivative-free coordinate search.
- The upstream `skewhypCheckPars` accepts `nu = 0` even though its error text says nu must be greater than zero and the density formulas are singular there. This port requires `nu > 0`.
- `ddskewhyp` is exposed as a stable central numerical derivative rather than carrying over the very large symbolic derivative expression from the R source.
