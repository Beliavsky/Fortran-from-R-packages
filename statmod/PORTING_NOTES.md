# Porting notes

## Scope

The goal is computational parity with `statmod` 1.5.2 while presenting a
Fortran-native numerical API. Plotting, printing, S3 dispatch, model-frame and
formula parsing, column/row names, and other R presentation machinery are not
reimplemented.

The translated numerical areas are:

- Digamma exponential-family functions;
- inverse-Gaussian density/CDF/quantile/RNG;
- Gaussian quadrature and probability-measure quadrature;
- expected unit-deviance means/variances, including the full upstream native
  Chebyshev coefficient tables;
- secure Gamma and negative-binomial GLM fitting;
- `fitNBP`;
- two-variance-component mixed models;
- normal and Gamma heteroscedastic REML scoring;
- GLM score testing;
- ELDA/limiting-dilution estimation and group/slope tests;
- matrix scaling, robust scale, forward selection, multiple/permutation tests,
  SAGE tests and growth-curve permutation tests;
- randomized quantile residual kernels for the supported families;
- Tweedie-family variance/link/deviance kernels.

## Use of r_mod

The supplied MIT-licensed `r_mod.f90` is reused for probability distributions,
RNGs, special functions and test helpers whenever applicable. New code was
added only where `statmod` needs behavior not supplied by that module, notably:

- the stable inverse-Gaussian tail formulas and inversion algorithm specific to
  `statmod`;
- the real-size negative-binomial summand needed by `expectedDeviance` (the
  supplied `r_mod` negative-binomial DPQ interface takes integer size);
- a small stable `expm1` helper used by ELDA because Fortran 2018 has no
  standard `expm1` intrinsic;
- LAPACK-backed weighted least-squares/SVD/SPD utilities needed by the model
  fitting algorithms.

No alternative implementation is added where an applicable `r_mod` routine is
already available.

## Native-code preservation

The upstream `src/gaussq2.f` algorithm was translated to free-form Fortran with
its Netlib/EISPACK attribution preserved. The full coefficient arrays and
piecewise Chebyshev logic from `src/expectedDeviance.c` were translated rather
than replaced by simulation or generic numerical integration.

## R-object boundaries

`mixedModel2`/`randomizedBlock`, `Digamma`, `tweedie`, `qresiduals`, and the
high-level ELDA routines use R formulas, family objects or S3 dispatch. Their
computational kernels are exposed directly in Fortran; the R object-system
shell is intentionally omitted. Integer group vectors replace R factors.

`remlscoregamma` is exposed for the upstream default log mean and log
dispersion links. The generic R link-object interface is not recreated.

## Optional Tweedie residual dependency

Upstream `qres.tweedie` calls the suggested R package `tweedie`. To retain this
computational path without an R runtime, this project vendors the previously
translated Tweedie computational core under `src/tweedie_dep/`. It retains its
own GPL-2.0-or-later notice. No second copy of `r_mod` is vendored.

## Validation

The included tests use the upstream package's saved regression output where it
provides stable deterministic references. In particular they check:

- inverse-Gaussian probabilities and quantiles;
- Legendre/Jacobi/normal quadrature nodes and weights;
- binomial, Gamma, negative-binomial and Poisson expected-deviance values;
- the published deterministic `mixedModel2Fit` example;
- a deterministic `glmnb.fit` example;
- ELDA equal-dose MLE identity;
- GLM score residualization, Holm adjustment, exact permutation/SAGE cases;
- `matvec`/`vecmat`, forward selection, robust scale, `fitNBP`, REML-Gamma and
  Tweedie residual smoke tests.

A clean direct GNU Fortran build uses Fortran 2018, runtime checking and
`-Werror=implicit-interface`. FPM is the supported project layout; if `fpm` is
not installed, the same source graph can be compiled directly and linked to
BLAS/LAPACK.
