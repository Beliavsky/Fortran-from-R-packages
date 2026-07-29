# evir-fortran

A self-contained modern Fortran translation of the computational algorithms in
R package **evir 1.7-4** for extreme-value analysis.

## Scope

The library includes:

- GEV and Gumbel block-maxima maximum-likelihood fitting;
- generalized Pareto ML and probability-weighted-moment fitting;
- nonhomogeneous Poisson point-process POT fitting;
- bivariate logistic POT fitting with fixed or jointly fitted margins;
- GEV return levels and profile-likelihood intervals;
- GPD quantiles, expected shortfall, Wald intervals, and profile intervals;
- GEV/GPD density, CDF, quantile, and random generation;
- declustering and threshold selection;
- Hill, mean-excess, extremal-index, shape-stability, quantile-stability,
  empirical-tail, QQ, record-development, and tail-curve diagnostics.

The original interactive plotting methods are represented by procedures that
return the numerical coordinates needed to make the plots. R S3 classes,
menus, graphics, `ts` attributes, and formatted printing are not reproduced.

## Build with FPM

```text
fpm build
fpm test
fpm run
fpm run --example block_maxima_example
```

A GNU Fortran build script is also supplied:

```text
sh scripts/build_and_test.sh debug
sh scripts/build_and_test.sh release
```

## Minimal example

```fortran
program example
    use evir, only : dp, gpd, riskmeasures, gpd_fit_result
    implicit none
    type(gpd_fit_result) :: fit
    real(dp) :: data(8), p(2), q(2), es(2)

    data = [0.2_dp, 0.5_dp, 0.9_dp, 1.2_dp, 1.7_dp, 2.4_dp, 3.1_dp, 5.0_dp]
    fit = gpd(data, threshold=1.0_dp, information='expected')
    p = [0.95_dp, 0.99_dp]
    call riskmeasures(fit, p, q, es)
    print *, q
    print *, es
end program example
```

## Dependencies

The numerical library has no external dependencies. It contains its own small
dense optimizer, numerical Hessian, matrix inversion, random-number generator,
and normal/chi-square quantile support.

## License

GPL-2.0-or-later. See `LICENSE`, `NOTICE.md`, and `licenses/`.
