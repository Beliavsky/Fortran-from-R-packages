# Porting notes

## Scope

This is a computational translation of R package `skewunit` 1.1 (dated
2026-07-31 in the supplied `DESCRIPTION`). All exported numerical/statistical
functionality is represented. `print.skewunit` is presentation-only and is
not translated as a library routine; ordinary Fortran I/O can print the
result types. The supplied package contains no plotting routines.

The original R sources and manual pages are retained under `upstream/` for
provenance and behavioral comparison.

## R-to-Fortran mapping

- `dasin/pasin/rasin` -> same names plus `rasin_vec`
- `dtriang/ptriang/rtriang` -> same names plus `rtriang_vec`
- `dUquad/pUquad/rUquad` -> `duquad/puquad/ruquad` plus `ruquad_vec`
- `dJSB/pJSB/rJSB` -> `djsb/pjsb/rjsb` plus `rjsb_vec`
- `dsbeta/psbeta/rsbeta` -> same names plus `rsbeta_vec`
- `dskewunit/pskewunit/rskewunit` -> same names plus vector helpers
- `estimate.skewunit` -> `estimate_skewunit`
- `choose.skewunit` -> `choose_skewunit`
- `cuberoot` -> `cuberoot`

Fortran is case-insensitive, so mixed-case R function names become lowercase
source identifiers naturally.

## Replacements for R dependencies

`stats::dbeta/pbeta/rbeta` are implemented with `log_gamma`, a continued
fraction for the regularized incomplete beta, and Gamma-ratio simulation.
`stats::dnorm/pnorm` are replaced by direct formulas and `erfc`.
`stats::optim` is replaced with standalone Brent and Nelder-Mead routines.
`stats::integrate` is replaced by specialized adaptive 15-point Gauss-Kronrod quadrature.
`pracma::hessian` is replaced by a central finite-difference Hessian.

The Hessian is at most 3 by 3, so its inverse is obtained with a small
pivoted Gauss-Jordan routine rather than an external BLAS/LAPACK dependency.

## Compatibility details

The shape-parameter routing of `dskewunit`, `pskewunit`, `rskewunit`, and the
estimator follows the intent and effective behavior of the R package:
`delta` is the sole shape when exactly one component is JSB/sbeta, while
`delta` and `delta2` are the two shapes when both components require one.

The R estimator deletes observations exactly equal to 0.5 when `family1` is
U-quadratic because that baseline density is zero there. The Fortran MLE
reproduces this behavior.

The upstream baseline density functions set their result to zero at and
outside the support even when `log=TRUE`. The translated baseline and skew
density functions preserve this unusual convention at the endpoints.

## Intentional numerical fixes/improvements

Three upstream implementation details are not copied literally:

1. In `pJSB`, `pUquad`, and `ptriang`, the R code applies `1-F` after taking a
   logarithm when both `lower.tail=FALSE` and `log.p=TRUE`. The Fortran API
   returns the conventional `log(1-F)`.
2. `rskewunit` in R accepts a proposal with probability `G/2`. The constant
   factor cancels in rejection sampling, but halves efficiency. The Fortran
   sampler accepts with probability `G`, producing the same target law.
3. The R rejection sampler contains a family-name typo in one shape dispatch
   condition. The Fortran version uses the intended `family2` shape routing.

The R code also performs two optimizer starts in several cases. The Fortran
port retains the two Nelder-Mead starts for multi-parameter models. Pure
one-dimensional cases use one deterministic bounded Brent search across the
full transformed interval, so a second start is unnecessary.

## CDF integration

A direct integral in `x` is numerically awkward when `family1=asin`, because
its density has integrable square-root singularities at both endpoints. The
port integrates after

```text
x = sin(theta)^2
```

which removes those singularities. This was necessary for high-accuracy
agreement with independent reference values.

## Random-number streams

The generated distributions match the upstream laws, but random sequences do
not match R's RNG stream. The port uses Fortran `random_number`, a polar normal
sampler, and Marsaglia-Tsang Gamma generation. `seed_skewunit_rng` provides a
repeatable seed within the Fortran implementation/compiler environment.

## Error handling

R uses `stop()` for invalid user arguments. Scalar Fortran distribution
routines generally return IEEE NaN for invalid shape/range parameters. The
estimator reports invalid input through `fit%convergence=2` and infinite
criterion values rather than terminating the process.

## FPM

The source layout is conventional FPM (`src/`, `test/`, `example/`) and the
manifest declares no external dependencies.
