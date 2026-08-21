# boot-fortran

Modern Fortran computational-core translation of the R `boot` package 1.3-32.

The upstream package implements the bootstrap algorithms from Davison and
Hinkley, *Bootstrap Methods and Their Application*.  This port keeps the
numerical/statistical algorithms and replaces R objects, formulas, S3 methods,
and graphics with array- and callback-based Fortran APIs.

## Implemented numerical areas

- ordinary, balanced, permutation, antithetic, importance, and balanced-importance resampling arrays
- conversion between bootstrap index and frequency arrays
- generic weighted nonparametric bootstrap driver
- normal, basic, percentile, studentized, BCa, and ABC confidence intervals
- normal-quantile interpolation used by upstream `boot.ci`
- infinitesimal, delete-one, positive, and regression empirical influence values
- linear bootstrap approximation
- exponential tilting
- importance weights, raw/ratio/regression moment estimators, tail probabilities, and quantiles
- weighted correlation, `var.linear`, `k3.linear`, `cum3`, logit, and inverse-logit helpers
- pointwise and overall bootstrap confidence envelopes
- frequency smoothing used by tilted bootstrap calculations
- fixed- and geometrically distributed block-bootstrap index generation
- empirical-likelihood and exponential-family empirical-likelihood numerical profiles
- likelihood-profile confidence limits
- general two-phase linear-programming solver corresponding to `simplex`
- simple multinomial saddlepoint density/CDF approximation
- censored case-resampling and product-limit/conditional-censor sampling primitives
- GLM diagnostic formulas from `glm.diag`
- nested-correlation bootstrap calculation

The facade module is `boot`.

## Example

```fortran
program bootstrap_mean
    use boot_kinds, only : dp
    use boot_core
    use boot_ci, only : percentile_ci
    implicit none

    real(dp) :: data(8,1), lo(1), hi(1)
    type(bootstrap_result) :: res

    data(:,1) = [2.0_dp, 3.0_dp, 5.0_dp, 7.0_dp, &
                 11.0_dp, 13.0_dp, 17.0_dp, 19.0_dp]

    call bootstrap_weighted(data, stat, 999, 'ordinary', res)
    call percentile_ci(res%t, [0.95_dp], lo, hi)

contains

    function stat(x, w) result(value)
        real(dp), intent(in) :: x(:,:), w(:)
        real(dp) :: value
        value = sum(x(:,1)*w)/sum(w)
    end function stat
end program bootstrap_mean
```

## Build

With FPM:

```text
fpm test
fpm run --example bootstrap_mean
```

The project has also been validated directly with GNU Fortran 2018 using
`-Wall -Wextra -Werror -fcheck=all`.

## License

The upstream `boot` package declares `License: Unlimited` and its source states
`Unlimited distribution is permitted`.  The complete upstream 1.3-32 snapshot
is retained under `upstream/boot-1.3-32`; this translation is distributed under
the same terms.

See `docs/TRANSLATION_COVERAGE.md` and `PORTING_NOTES.md` for the exact scope.
