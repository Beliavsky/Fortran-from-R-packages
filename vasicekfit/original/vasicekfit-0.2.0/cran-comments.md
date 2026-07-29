## Test environments

* local Linux (Rocky 9.7), R 4.5.3
* R-hub (Windows, R-devel) via GitHub Actions

## R CMD check results

0 errors | 0 warnings | 0 notes

## Submission notes

This is an update of an existing CRAN package (0.1.0 -> 0.2.0).

User-visible changes since 0.1.0 (see NEWS.md):

* `predict(type = "response")` now returns the conditional mean loss rate
  (the effective PD), matching the convention of `predict.glm()`. The previous
  behaviour (conditional median) is available via `alpha = 0.5`.
* `vcov()`, `confint()`, and `summary()` gain a `type` argument offering HAC
  standard errors via the suggested `sandwich` package; the default is
  unchanged.
