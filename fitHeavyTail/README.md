# fitHeavyTail for modern Fortran

A self-contained modern Fortran translation of the computational routines in
**fitHeavyTail 0.2.0.9000**, packaged for the Fortran Package Manager (FPM).
The original R package estimates locations, scatter matrices, covariance
matrices, skewness vectors, and Student-t degrees of freedom under heavy tails.

## Features

- Tyler's shape estimator with spatial-median initialization and scale recovery
- Multivariate Cauchy maximum-likelihood estimation
- Multivariate Student-t EM/PX-EM estimation
- Fixed, one-shot, or iterative degrees-of-freedom estimation
- OPP, harmonic OPP, and the documented POP variants
- ECM and ECME degrees-of-freedom updates
- Observation weights
- Low-rank factor-analysis covariance constraints
- IEEE-NaN row removal or conditional-moment EM imputation
- Optional finite-sample covariance scaling
- Generalized-hyperbolic multivariate skewed-t EM/PX-EM estimation
- Real-order modified Bessel K evaluation and ratios
- Marginal-kurtosis, cross-cumulant, Hill, and Pareto tail estimators
- No external FPM or numerical-library dependencies

Plotting, R object infrastructure, vignettes, and saved R regression fixtures
are intentionally omitted. The original computational R sources, manual pages,
metadata, and license are retained under `original/`.

## Build with FPM

```text
fpm build
fpm test
fpm run
fpm run --example tyler_cauchy_example
fpm run --example factor_missing_example
fpm run --example skew_t_example
```

The package requires a Fortran 2018 compiler.

## Minimal example

```fortran
program robust_t_example
   use fitheavytail
   implicit none

   real(dp) :: x(300,3)
   type(heavy_tail_fit) :: fitted

   call random_mvt_identity(300,3,6.0_dp,x,12345)
   call fit_mvt(x,fitted,nu_method='iterative', &
      nu_iterative_method='POP')

   print '(a,f10.4)', 'estimated nu: ',fitted%nu
   print '(*(f12.6,1x))',fitted%covariance
end program robust_t_example
```

## Main API

| R name | Fortran procedure |
|---|---|
| `fit_Tyler` | `fit_tyler` |
| `fit_Cauchy` | `fit_cauchy` |
| `fit_mvt` | `fit_mvt` |
| `fit_mvst` | `fit_mvst` |
| `nu_OPP_estimator` | `nu_opp_estimator` |
| `nu_POP_estimator` | `nu_pop_estimator` |

R strings such as `"POP-approx-2"`, `"ECM"`, `"MLE-diag"`, and
`"cross-cumulants"` are retained where they select numerical methods.

Fortran has no universal missing-value scalar. Missing observations are
represented with IEEE quiet NaNs. Pass `na_rm=.false.` to `fit_mvt` to use the
conditional-moment EM calculations rather than dropping incomplete rows.

See [API.md](API.md), [PORTING.md](PORTING.md), and
[TRANSLATION_COVERAGE.md](TRANSLATION_COVERAGE.md) for details.

## Testing

`fpm test` builds four numerical test programs. Strict GNU Fortran scripts are
also included:

```text
./run_gfortran_tests.sh
run_gfortran_tests.bat
```

## License

The original package declares GPL-3. This translation is distributed under
**GPL-3.0-only**. See [LICENSE](LICENSE) and [NOTICE.md](NOTICE.md).
