# Porting notes

## Numerical model

The port preserves the eight `poweRlaw` probability models and truncation at
`xmin`.  Discrete data are stored as `real(dp)` for a uniform API but are
validated to be integer-valued when attached to a discrete model.

For the discrete power law, the normalizer is the Hurwitz zeta
`zeta(alpha,xmin)`.  Calling the supplied `pracma::zeta` implementation during
thousands of likelihood evaluations is very slow near alpha=1 because it uses
a direct eta-series.  The Fortran port therefore evaluates Hurwitz zeta by
Euler-Maclaurin summation in that regime and uses the supplied pracma zeta in
the rapidly convergent high-alpha, xmin=1 case.  This also avoids subtracting
large nearly equal zeta/partial-sum values when xmin is large.

## Estimation

The continuous power-law, continuous exponential, and discrete exponential
MLEs have closed forms, so the Fortran port evaluates those formulas directly.
The R package routes them through `optim`; the optimum is the same apart from
optimizer tolerance.

Discrete power-law and truncated-Poisson one-parameter fits use bounded golden
section search on transformed parameters.  Lognormal and Weibull fits use a
standalone two-dimensional Nelder-Mead optimizer.  Candidate parameter grids
are supported by `estimate_pars` and `estimate_xmin`.

## xmin and goodness of fit

`estimate_xmin` implements the upstream Clauset-style scan over candidate lower
cutoffs, refitting parameters at each candidate and minimizing either the KS
or reweighted distance.  The discrete empirical CDF is evaluated on the
integer support; the continuous empirical CDF uses `(0:n-1)/n`, matching the
upstream implementation.

## Tail probabilities

The R sources contain several family-specific upper-tail conventions and
normalizers that are not mutually consistent (notably the custom discrete
power-law upper tail and some truncated discrete wrappers).  The Fortran
object API uses a consistent conditional upper tail, `1-CDF`, while preserving
upstream lower-tail behavior used by fitting and KS estimation.  This is a
small compatibility correction rather than a change to the fitted likelihoods.

## Bootstrap

`bootstrap` implements the ordinary nonparametric resampling workflow.
`bootstrap_p` implements the semi-parametric lower-body plus fitted-tail
workflow.  `xmins`, parameter grids, `xmax`, distance choice, and seed controls
are supported.  Parallel R orchestration and timing/progress messages are
omitted; v0.1.0 runs simulations serially.

The intrinsic Fortran RNG is used, so a given integer seed is reproducible
within this port but does not reproduce R's RNG stream.

## R infrastructure omitted

Reference-class/S4 machinery, plotting, printing/show methods, and the
parallel-cluster layer are not computational algorithms and are omitted.
