# sadists-fortran

Modern Fortran translation of the computational code in Steven E. Pav's R
package `sadists` 0.2.6.

Version 0.2.0 refactors the approximation layer to use the standalone
`pdqutils-fortran` package instead of maintaining a private copy of PDQutils
algorithms.  The public sadists API remains compatible with v0.1.0.

The library implements the package's 12 probability-distribution families and
all 48 density/CDF/quantile/random-generation operations.  The R/Shiny
`runExample()` application is intentionally not translated.

## Requirements

- Fortran 2018 compiler
- FPM (Fortran Package Manager)
- `pdqutils-fortran` 0.1.0 (vendored under `vendor/` for a self-contained build)

There are no runtime dependencies on R, hypergeo, or orthopolynom.

## Build and test

```text
fpm build
fpm test
fpm run --example basic
```

The FPM manifest uses the local dependency:

```toml
[dependencies]
pdqutils-fortran = { path = "vendor/pdqutils-fortran" }
```

## Distribution families

- doubly noncentral beta
- doubly noncentral eta
- doubly noncentral F
- doubly noncentral t
- K-prime
- lambda-prime
- upsilon
- weighted sums of powers of noncentral chi-squares
- weighted sums of log noncentral chi-squares
- products of powers of noncentral chi-squares
- products of doubly noncentral F variates
- products of normal variates

The public Fortran names follow the R package: `ddnf`, `pdnf`, `qdnf`,
`rdnf`, etc. Random generators are subroutines whose first argument is the
output array. For example:

```fortran
program demo
    use sadists
    implicit none
    real(dp) :: x(1000)

    call rdnf(x, 40.0_dp, 80.0_dp, 1.5_dp, 2.5_dp)
    print *, sum(x)/real(size(x),dp)
    print *, pdnf(1.0_dp, 40.0_dp, 80.0_dp, 1.5_dp, 2.5_dp)
end program demo
```

For multi-component distributions, parameters such as weights, degrees of
freedom, noncentralities, powers, normal means, and standard deviations are
rank-one arrays. Different parameter-array lengths are recycled in the same
way as the upstream R `mapply` calls.

R argument names map as follows:

- `log` -> `log_density`
- `lower.tail` -> `lower_tail`
- `log.p` -> `log_p`
- `order.max` -> `order_max`

The four doubly-noncentral scalar-parameter families also provide `_vec`
helpers for array evaluation of density, CDF, and quantile values.

## Numerical architecture

`sadists` computes most approximate density/CDF/quantile functions from raw
cumulants. In v0.2.0 those cumulants are passed to the canonical standalone
PDQutils-fortran implementation:

```text
sadists parameters
      |
      v
sadists moments/cumulants
      |
      v
pdqutils-fortran
 Edgeworth / Cornish-Fisher
      |
      v
density / CDF / quantile
```

`sadists_approximations.f90` is now only a compatibility adapter mapping the
v0.1 internal names to `dapx_edgeworth`, `papx_edgeworth`, `qapx_cf`, and
`as269`. `sadists_special.f90` similarly delegates common normal probability
utilities and moment/cumulant conversion to PDQutils while retaining
sadists-specific chi-square moments, special functions, and simulation code.

These remain approximations, just as in the R package. In particular, an
Edgeworth CDF need not be globally monotone for difficult parameter choices or
aggressive expansion order. See `PORTING_NOTES.md`.

## Source layout

```text
src/sadists_kinds.f90           precision/constants
src/sadists_special.f90         sadists-specific special functions, moments, RNGs
src/sadists_approximations.f90  thin PDQutils compatibility adapter
src/sadists_distributions.f90   translated distribution families
src/sadists.f90                 public umbrella module
vendor/pdqutils-fortran/        standalone PDQutils-fortran dependency
test/                           regression, RNG, and dependency-integration tests
example/                        FPM example
upstream/                       retained upstream sadists/provenance sources
```

## Licensing

Both the sadists-derived translation and the vendored PDQutils-fortran
dependency are LGPL-3.0-or-later. See `LICENSES.md`, `NOTICE.md`, and the
license/provenance files under `vendor/pdqutils-fortran/` and `upstream/`.
