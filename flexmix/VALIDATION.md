# Validation

Validation was performed on Linux with GNU Fortran 14.2.0. The package uses only
standard Fortran/FPM facilities and has no external numerical-library links, avoiding
separately installed BLAS/LAPACK requirements on Windows.

## Compiler validation

All maintained package modules were compiled together with warnings treated as errors in
both configurations:

```text
-O0 -g -std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -fbacktrace
-O2 -std=f2018 -Wall -Wextra -Wpedantic -Werror
```

In both configurations:

- `test/test_flexmix.f90` completed with `All flexmix tests passed.`
- `example/basic_flexmix.f90` ran successfully and recovered the intended two-component
  deterministic Gaussian regression mixture.
- `validation/parity_driver.f90` and `validation/independent_parity.py` completed
  successfully.

The deterministic suite covers the translated EM engine; Gaussian, Poisson, binomial,
Gamma inverse-link, multinomial, conditional-logit, robust-background, structural-zero,
fixed/shared-coefficient, elastic-net, mixed-effects, censored fixed/mixed-effects, and penalized-smooth regression
models; numeric offsets; multivariate normal, factor-analyzer, binary, Poisson, and mixed
binary/Gaussian mixtures; positive univariate mixtures; concomitant priors; grouped
observations; accessors; effective-df information criteria; relabeling and component
removal; parameter packing/replacement; analytic scores; observed-information
covariance/standard errors; simulation; candidate-k model selection; regression bootstrap;
bootstrap LR testing; and deterministic example-data generators.

## Independent numerical parity

`validation/parity_driver.f90` was compiled with the runtime-check configuration and
`validation/independent_parity.py` independently recomputed reference results with NumPy,
SciPy, scikit-learn, and statsmodels. The checks passed for:

- Gaussian weighted least-squares coefficients and residual scale;
- Poisson and binomial IRLS coefficients;
- multivariate-normal means and unbiased covariance;
- Gamma inverse-link coefficients and FlexMix's deviance-based shape estimate;
- Poisson regression with a numeric offset;
- multinomial intercept-only log-odds/class probabilities;
- fixed/shared Gaussian-regression coefficients;
- one-choice-per-stratum conditional-logit coefficients using an independent SciPy
  likelihood optimizer;
- factor-analysis uniquenesses using an independent bounded SciPy optimization of the
  upstream concentrated maximum-likelihood objective;
- fixed-lambda Gaussian elastic-net/lasso coefficients against scikit-learn;
- grouped random-intercept mixed-model fixed effects, variance components, and marginal
  log likelihood against statsmodels MixedLM;
- left-censored fixed-effect Gaussian regression parameters and observed-data log likelihood
  against a separately optimized SciPy censored-normal likelihood;
- left-censored random-intercept mixed-model parameters and observed-data log likelihood
  against a separately optimized SciPy likelihood, plus a deterministic bivariate-normal
  integration check for the multiple-censor QMC path; and
- fixed-smoothing-parameter Gaussian, Poisson, and binomial smooth regression coefficients,
  with the Gaussian path checked against the closed-form penalized solve and the GLM paths
  checked against independent SciPy penalized-likelihood optimizers.

The generated comparison values are retained in `validation/fortran_parity.csv`.

## Static release audit

`python validation/audit_source.py` reports:

```text
SOURCE AUDIT PASSED: 17 Fortran files; 66/68 coverage mapping consistent.
```

The audit checks explicit dummy `INTENT`/`VALUE` attributes, one dummy declaration per
line with trailing FORD `!!` documentation, semicolon-separated statements, legacy real
kinds/D-exponents, self-comparison NaN tests, unsafe fast-math flags, duplicate Fortran
sources, build/archive products, vendored dependency directories, coverage arithmetic,
mapping/source existence, and README/manifest agreement. All maintained Fortran source
also fits gfortran's default free-form line length without `-ffree-line-length-none`.

## FPM availability limitation

FPM is not installed in the execution environment used for this translation. Immediately
before packaging, the required commands were explicitly attempted and each returned shell
status 127 (`fpm: not found`):

```text
fpm build
fpm test
fpm clean --all
```

Therefore this document does **not** claim those three FPM commands ran successfully here.
Instead, the same source set, test main, and example main declared in `fpm.toml` were built
and run directly with gfortran in clean external directories, in both debug/runtime-check
and optimized configurations. The TOML manifest was parsed successfully with Python's
`tomllib`.

No compiler output is created inside the package tree by this validation workflow, so the
final source tree is manually clean even though `fpm clean --all` cannot run on this host.

## Upstream archive

The supplied upstream `flexmix-master.zip` used for this translation has SHA-256:

```text
55539189569e874fd78652edc308e7a8b71cd194f72814ca832a591509a29d20
```
