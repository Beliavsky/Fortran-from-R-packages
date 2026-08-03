# segmented-fortran

A self-contained modern Fortran/FPM implementation of the computational core of
R's `segmented` package, including its `nlme`-based mixed-effects path.

## Implemented model families

- Continuous segmented (broken-line) Gaussian regression.
- Piecewise-constant stepmented Gaussian regression.
- Binomial-logit and Poisson-log segmented and stepmented regression.
- Segmented linear mixed-effects models using the translated `nlme` covariance,
  likelihood, and BLUP implementation.
- Multiple estimated breakpoints, including multiple breakpoints on one
  covariate by repeating its column in `z`.
- Observation weights and offsets for ordinary and generalized models.

Continuous models use the package's hinge/indicator linearization. Step models
refit thresholds over intervals between observed covariate values. Breakpoint
updates are clamped to configurable quantile limits and safeguarded by line
search.

## Supporting computations

- Segmented/stepmented model matrices and prediction.
- Segment slopes and intercepts.
- Delta-style breakpoint standard errors and normal confidence intervals.
- Average annual percent change/average slope (`aapc`).
- Broken-line evaluation.
- Grid Davies-style and fixed-breakpoint score tests.
- Fixed-breakpoint power approximation.
- BIC selection over zero through a requested maximum number of breakpoints.
- Numeric-vector convenience wrappers and R-recognizable aliases such as
  `segmented_lm`, `segmented_glm`, `segmented_lme`, `stepmented_lm`, `segreg`,
  and `stepreg`.

## Basic use

```fortran
use segmented

type(segmented_result) :: fit
real(dp) :: psi0(1)

psi0 = 5.0_dp
call segmented_lm(y, x, z, psi0, fit)
print *, fit%breakpoints
print *, fit%coefficients
```

For mixed effects:

```fortran
use segmented

type(segmented_lme_result) :: fit
type(segmented_lme_options) :: options

call segmented_lme(y, x, z, psi0, random_design, group, fit, options)
```

## Build

```text
fpm build
fpm test
fpm run demo_segmented
```

Direct GNU Fortran scripts are included for environments without FPM:

```text
./run_gfortran_tests.sh debug
./run_gfortran_tests.sh release
```

## Scope differences from R

This is a numerical library, not a reproduction of R's formula and S3 runtime.
Formula parsing, plotting, printing, interactive history drawing, `lm`/`glm`
object mutation, and bundled `.rda` datasets are not part of the compiled API.
The bootstrap-restart, constrained-slope, ARIMA, and full random-breakpoint
parameterizations are not reproduced exactly. The mixed-effects routine
supports explicit fixed and random design matrices and uses the translated
`nlme` engine. See `PORTING.md` and `TRANSLATION_COVERAGE.md`.
