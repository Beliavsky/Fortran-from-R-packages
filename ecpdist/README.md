# ecpdist-fortran

Modern Fortran 2018 / FPM translation of the computational code in the R
package **ecpdist 0.2.1**, which implements the Extended Chen-Poisson (ECP)
lifetime distribution.

Upstream authors: Ana Abreu, Ivo Sousa-Ferreira, and Cristina Rocha.

Reference:

> Sousa-Ferreira, I., Abreu, A.M. & Rocha, C. (2023). The Extended
> Chen-Poisson Lifetime Distribution. REVSTAT - Statistical Journal 21(2),
> 173-196. DOI: 10.57805/revstat.v21i2.405.

## Implemented numerical API

The top-level module is `ecpdist`.

- `decp` - density or log-density
- `pecp` - CDF/survival probability, with lower/upper tail and log probability
- `qecp` - quantile, with lower/upper tail and log-probability input
- `recp` - pseudo-random sample generation by inversion
- `secp` - survival/CDF and cumulative-hazard compatibility interface
- `hecp` - hazard or log-hazard
- `ecp_cumhaz` - cumulative hazard
- `ecp_kmoment` - k-th raw moment with numerical error estimate
- `ecp_kmoment_cond` - conditional k-th raw moment given `X > x`
- `ecp_mrl` - mean residual life
- `ecp_shape` - Bowley skewness or Moors kurtosis

The R plotting helper `ecp_plot` is intentionally omitted. The quantities it
plots are all available from the numerical functions above.

## Build

```sh
fpm build
fpm test
fpm run --example basic
```

A compiler-only validation can be performed with GNU Fortran using strict
Fortran 2018 checks. See `PORTING_NOTES.md` for the command used for this
release.

## Example

```fortran
program demo
    use ecpdist, only: dp, decp, pecp, qecp
    implicit none

    print *, decp(1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
    print *, pecp(1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
    print *, qecp(0.5_dp, 1.0_dp, 1.0_dp, 1.0_dp)
end program demo
```

## Numerical implementation notes

Probability calculations use log-survival algebra so both positive and
negative `phi` remain stable without forming avoidable large exponentials.
Raw and conditional moments are evaluated by adaptive 15-point Gauss-Kronrod
integration over the quantile representation. This computes the same
expectations as the upstream `stats::integrate()` formulas while avoiding their
endpoint singularity at `y=0`.

See `PORTING_NOTES.md` for two corrected upstream `qecp` option bugs.

## License

The upstream package declares `License: GPL-3`. This translation is therefore
distributed under **GPL-3.0-only**. See `LICENSE` and `LICENSES.md`.
