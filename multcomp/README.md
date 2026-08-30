# multcomp-fortran

Modern Fortran translation of the computational core of R `multcomp` 1.4-32.
The project provides general linear-hypothesis tests, simultaneous confidence
intervals, multiple-testing adjustments, contrast construction, compact-letter
displays, and multiple-marginal-model covariance assembly without requiring R.

## Main API

```fortran
use multcomp
```

The public numerical API includes:

- `make_parm` -- construct a coefficient/covariance/df parameter object.
- `parm_from_iid` -- construct covariance information from IID influence values.
- `contr_mat` -- all ten upstream `contrMat()` families: Dunnett, Tukey,
  Sequen, AVE, Changepoint, Williams, Marcus, McDermott, UmbrellaWilliams,
  and GrandMean.
- `glht_fit` -- test arbitrary hypotheses `K beta = rhs`.
- `glht_identity` and `glht_coefficients` -- all- or selected-coefficient tests,
  covering the computational role of `cftest()`.
- `glht_test` -- univariate, single-step, free, Shaffer, Westfall, Bonferroni,
  Holm, Hochberg, Hommel, BH/FDR, BY, or no p-value adjustment.
- `glht_confint` -- univariate or simultaneous normal/t confidence intervals.
- `glht_global_test` -- generalized-inverse Wald chi-square or F tests.
- `glht_critical_value` -- simultaneous or univariate critical constants.
- `compact_letter_display`, `cld_from_tukey_test`, and
  `cld_from_tukey_confint` -- Piepho insert-absorb compact-letter displays.
- `mmm_parm_from_iid` -- concatenate marginal models and recover cross-model
  covariance from IID estimating-function contributions.
- `block_diagonal_matrix` -- native equivalent of the block-diagonal utility
  used by `mmm`/`mlf`.

`probability_control`, `genz_bretz`, `tvpack`, and `miwa` from the vendored
`mvtnorm-fortran` backend are re-exported through the public `multcomp` module.

## Example

```fortran
use multcomp

type(parm_type) :: parameters
type(contrast_matrix_type) :: tukey
type(glht_type) :: hypotheses
type(mtest_type) :: tests
real(dp) :: beta(3), covariance(3, 3)

beta = [10.0_dp, 12.0_dp, 15.0_dp]
covariance = 0.0_dp
covariance(1, 1) = 1.0_dp
covariance(2, 2) = 1.0_dp
covariance(3, 3) = 1.0_dp

call make_parm(beta, covariance, parameters, df=30.0_dp)
call contr_mat([12.0_dp, 12.0_dp, 12.0_dp], 'Tukey', tukey)
call glht_fit(parameters, tukey%value, hypotheses)
call glht_test(hypotheses, 'single-step', tests)
```

See `example/simultaneous_inference.f90` for a complete program.

## Validation

`tools/run_strict_tests.sh` builds the vendored `mvtnorm` backend, the full
library, all regression tests, and the public example with GNU Fortran using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fimplicit-none
-fcheck=all -ffpe-trap=invalid,zero,overflow
```

The tests include direct numerical fixtures from the upstream saved R outputs:

- the one-sided `parm`/`glht` simultaneous-confidence-interval example;
- the cholesterol all-pair Tukey example for univariate, Shaffer, and Westfall
  p-values;
- the warpbreaks Tukey single-step p-values and simultaneous critical value;
- published `contrMat()` matrices;
- all R `p.adjust` families used by `multcomp`;
- compact-letter and multiple-marginal-model covariance invariants.

## Native-language adaptations

The statistical inference engine is translated. R formula parsing, S3 methods,
model-frame/model-matrix extraction, model-specific `coef()`/`vcov()` methods,
and plotting/printing are intentionally not reproduced. A native caller passes
coefficients, covariance, degrees of freedom, and contrast matrices directly.

`mcp()` is therefore represented by explicit contrast matrices plus
`contr_mat()`, rather than by parsing R model terms. For `mmm`, the caller
supplies IID influence values directly rather than relying on `sandwich::estfun`.
See `TRANSLATION_NOTES.md` for detailed parity qualifications.

## Build with FPM

```text
fpm build
fpm test
fpm run --example simultaneous_inference
```

The project is GPL-2.0-only, matching upstream `multcomp` and the vendored
`mvtnorm` dependency.
