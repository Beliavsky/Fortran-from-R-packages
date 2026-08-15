# Validation

Compiler: GNU Fortran 14.2.0

Core validation flags:

```text
-std=f2018 -Werror=implicit-interface -Werror=trampolines -fcheck=all -O0
```

Fresh-source validation compiles:

1. supplied `splines-fortran` source;
2. supplied `nlme` source;
3. supplied `survival-fortran` source;
4. vendored `gamlss.dist-fortran v0.3.0` source;
5. all `gamlss-fortran` modules;
6. ten integration tests and nine examples.

Expected test output:

```text
test_core: PASS
test_lms_selection: PASS
test_v02: PASS
test_v03: PASS
test_v04: PASS
test_v05: PASS
test_v06: PASS
test_v07: PASS
test_v08: PASS
test_v09: PASS
```

`test_v02` exercises:

- doubly censored normal regression and recovery of location/scale;
- grouped random intercepts with the supplied `nlme` implementation used for
  the initial Gaussian variance ratio;
- fractional-polynomial search, including the repeated-power rule;
- 1D LOESS reconstruction and effective degrees of freedom;
- varying-coefficient P-spline construction;
- monotone P-spline inequality enforcement;
- adaptive Lp category fusion;
- forward GAIC/BIC column selection;
- coefficient profile likelihood.

The basic example recovers the synthetic normal location model approximately as

```text
Converged: T
mu coefficients:     2.0005    1.4103
sigma:                0.2477
global deviance:      5.5971
```

The v0.2 example reports:

```text
Selected fractional-polynomial power:   -1.00
Weighted residual deviance:   0.016105
```

All project/test/example Fortran lines are kept within the standard 132-column
free-form limit.


`test_v03` exercises:

- correlated random intercept/slope covariance estimation and fixed-effect recovery;
- delayed-entry exact likelihoods plus interval2/counting-process adapters;
- additive persistent P-spline composition and 2D tensor penalties;
- backward and bidirectional BIC/GAIC column selection;
- case-resampling bootstrap percentile intervals;
- profile-likelihood confidence intervals.

The v0.3 example reports approximately:

```text
Random-slope fixed coefficients:     1.0291    0.7500
Random covariance (v11,v12,v22):    0.01917   0.00703   0.00687
Delayed-entry normal (mu,sigma):     0.4271    1.0824
```
`test_v04` exercises:

- exact Gaussian/NO AR(1) residual-correlation fitting through the supplied
  `nlme` GLS backend;
- simultaneous random intercepts on `mu` and `sigma`;
- bidirectional BIC selection on the `sigma` equation;
- worm-plot cubic coefficients, influence measures and Jarque-Bera diagnostics.

The v0.4 example reports approximately:

```text
Correlated-NO coefficients:     0.8024    1.1544
Fitted AR(1) correlation:      0.4326
```

The clean v0.4 link emits no GNU executable-stack/trampoline warnings.



## v0.5.0

`test_v05` checks:

- Gamma GAMLSS fitting with a fixed AR(1) working correlation;
- simultaneous random intercept/slope blocks on `mu` and `sigma`;
- three-fold matrix-first held-out log-score comparison;
- full-family randomized quantile residuals using SHASH.

`example/v05_extended.f90` simulates a correlated Gamma response and reports
correlated-RS coefficient recovery plus a three-fold held-out mean log score.


## v0.6.0

`test_v06` checks both limitations called out in v0.5:

- a Gamma response generated from an AR(1) Gaussian copula is fitted with the
  exact continuous-margin copula joint likelihood; the marginal slope and AR(1)
  dependence are recovered and the copula likelihood contribution is positive;
- a normal GAMLSS with random effects on both `mu` and `sigma` estimates a
  positive cross-parameter covariance rather than separate block covariances.

A separate rank-2 smoke test uses random intercept/slope designs on both `mu`
and `sigma`; the fitter returns the expected 4 by 4 joint covariance, including
cross-parameter intercept and slope covariance entries.

`example/v06_extended.f90` reports approximately:

```text
Copula Gamma mu coefficients:     0.4767    0.5953
Gaussian-copula AR(1) rho:        0.6923
Copula log-likelihood contribution: 23.9191
Joint-RE fixed mu coefficients:   1.0732    0.6175
Estimated Cov(b_mu,b_sigma):      0.03294
```

The generating copula AR(1) correlation is 0.55 in the example.


## v0.7.0

`test_v07` checks:

- the MVN rectangle integrator against the exact bivariate quadrant identity
  `P(X<=0,Y<=0)=1/4+asin(rho)/(2*pi)`;
- endpoint atom/left-limit handling for the mixed-mass `BEINF` family;
- grouped NBI data generated from an AR(1) Gaussian copula and fitted through
  the discrete rectangle likelihood;
- a two-dimensional joint `mu`/`sigma` random-effect model fitted by direct
  five-point Gauss-Hermite marginalization, including an SPD full covariance.

`example/v07_extended.f90` reports approximately:

```text
Discrete-copula fitted AR(1) rho:    0.4794
GHQ marginal log likelihood:       -13.6400
GHQ covariance (v11,v12,v22):       0.00357  -0.00097   0.00026
```

The generating discrete-copula AR(1) correlation is 0.50.


`test_v08` exercises:

- AIS against the closed-form Gaussian random-intercept marginal likelihood;
- a six-dimensional joint `mu`/`sigma` random-effect block, beyond the v0.7
  tensor-GHQ dimension cap;
- finite group posterior covariance and effective-sample-size diagnostics;
- optional outer AIS parameter refinement without recursive optimizer calls.

The v0.8 example reports approximately:

```text
AIS latent dimension: 6
AIS marginal log likelihood:       9.4983
Minimum group ESS:      259.3
```


`test_v09` exercises the upstream `getMarginal()` translation.  For a log-link
parameter with Gaussian random-effect standard deviation `s`, the integration
method is checked against the analytic identity
`E[exp(eta + U)] = exp(eta + s**2/2)`.  The suite also verifies identity-link
marginalization, the 999-point quantile rule, Monte Carlo averaging,
parameter-specific inverse links, the random-intercept object adapter, and the
stored random-effect EDF/`sigma_b` summary.

The v0.9 example reports approximately:

```text
Conditional means:              1.0000    2.0000    4.0000
Marginal means (integrate):     1.1972    2.3944    4.7889
Marginal means (qfunction):     1.1940    2.3880    4.7760
```
