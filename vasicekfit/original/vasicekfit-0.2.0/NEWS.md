# vasicekfit 0.2.0

* `predict(type = "response")` with `alpha = NULL` now returns the conditional
  mean loss rate (the effective PD), `Phi(qnorm(p) + sum kappa_j u_j)`, matching
  the convention of `predict.glm()`. Previously it returned `pnorm()` of the
  linear predictor, which is the conditional *median*; that quantity is now
  available via `alpha = 0.5`.
* `vcov()`, `confint()`, and `summary()` gain a `type` argument. `type = "HAC"`
  computes heteroskedasticity- and autocorrelation-consistent standard errors
  via `sandwich::lrvar()`; extra arguments are forwarded to it. The default
  `type = "iid"` is unchanged. `sandwich` is a suggested dependency.

# vasicekfit 0.1.0

* Initial CRAN release.
