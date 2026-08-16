# pdqutils-fortran

Modern Fortran 2018 translation of the computational code in Steven E. Pav's
R package **PDQutils 0.1.6**.

The library implements moment/cumulant conversion and probability density,
distribution, quantile, and random-generation approximations based on
Edgeworth, generalized Gram-Charlier, and Cornish-Fisher expansions.  It is a
standalone FPM package and has no R runtime dependency.

## Implemented upstream exports

- `AS269` -> `as269`, `as269_orders`, `as269_vector`
- `moment2cumulant`
- `cumulant2moment`
- `dapx_edgeworth`
- `papx_edgeworth`
- `dapx_gca`
- `papx_gca`
- `qapx_cf`
- `rapx_cf`

Array helpers with `_vec` suffix are also supplied for the scalar PDQ
functions.

## Generalized Gram-Charlier bases

`dapx_gca` and `papx_gca` support all bases exposed by the R package:

```fortran
use pdqutils, only : gca_normal, gca_gamma, gca_beta, &
                    gca_arcsine, gca_wigner
```

The orthogonal polynomials are evaluated from standalone Hermite,
generalized-Laguerre, and Jacobi recurrences.  Consequently the Fortran
package does not require the R `orthopolynom` or `moments` packages, nor a
Fortran translation of either package.

For gamma bases, explicit parent parameters can be supplied with `shape=` and
`scale=`.  For beta bases use `shape1=` and `shape2=`.  When omitted, the
parameters are inferred from the first two moments, as in PDQutils.

## Example

```fortran
program demo
    use pdqutils, only : dp, dapx_gca, papx_gca, gca_gamma
    implicit none
    real(dp) :: m(6)
    integer :: k, j

    ! Moments of chi-square(30).
    do k=1,6
        m(k)=1.0_dp
        do j=0,k-1
            m(k)=m(k)*2.0_dp*(15.0_dp+real(j,dp))
        end do
    end do

    print *, dapx_gca(30.0_dp,m,basis=gca_gamma,support_lo=0.0_dp)
    print *, papx_gca(30.0_dp,m,basis=gca_gamma,support_lo=0.0_dp)
end program demo
```

Build and test with FPM:

```text
fpm test
fpm run --example basic
```

## Numerical caveat

These routines reproduce approximation methods, not exact probability laws.
As documented by upstream PDQutils, Edgeworth and Gram-Charlier CDF
approximations are not guaranteed to be monotone, and Cornish-Fisher
quantiles are not guaranteed to be monotone.  The density and CDF routines
follow upstream behavior by clipping negative densities to zero and CDF
values to `[0,1]`.

See `PORTING_NOTES.md` for implementation details and compatibility fixes.

## License

PDQutils is LGPL version 3 or later according to its source headers.  The
translated code retains that licensing.  See `LICENSE`, `LICENSES.md`, and
`NOTICE.md`.  The complete supplied R source is retained under `upstream/` for
provenance.

## Validation

The release was compiled with GNU Fortran 14.2 using both an optimized build
and a runtime-checked Fortran 2018 build with implicit-interface diagnostics.
The shipped tests report:

```text
test_edgeworth:   PASS
test_gca:         PASS
test_moments_cf:  PASS
test_rng_vector:  PASS
```

FPM itself was not installed in the translation environment, so the same
`src/`, `test/`, and `example/` layout was compiled directly with gfortran.
The TOML manifest was parsed independently before packaging.
