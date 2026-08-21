# dirmult-fortran

Modern Fortran/FPM translation of the computational code in the R package
`dirmult` 0.1.3-5 by Torben Tvedebrink.

The upstream package estimates Dirichlet-multinomial parameters, computes
observed/expected Fisher information, evaluates profile likelihoods, fits a
common overdispersion parameter across multiple count tables, and provides
simulation utilities. This port translates those numerical routines and omits
R-specific object/data-frame presentation machinery.

## License and provenance

The upstream `DESCRIPTION` declares `License: GPL (>= 2)`. This translation is
therefore distributed under GPL-2.0-or-later. The GPL v2 text is in `LICENSE`.
The exact user-supplied upstream archive is preserved as
`provenance/dirmult-master.zip`, and selected upstream source/metadata is in
`provenance/upstream/`.

## Build

With Fortran Package Manager:

```text
fpm build
fpm test
fpm run --example basic_fit
```

The source requires Fortran 2018 only for the build/test mode used during
validation; the code itself uses broadly supported modern Fortran features.
There are no external numerical-library dependencies.

## Main API

Use the umbrella module:

```fortran
use dirmult_fortran
```

Important procedures are:

- `fit_dirmult`: Dirichlet-multinomial MLE by Fisher scoring, with observed or
  expected Fisher information.
- `summarize_dirmult`: MLE/MoM estimates and standard errors.
- `weir_mom`: Weir-Hill method-of-moments estimate of theta and its standard
  error.
- `dirmult_loglik`: parameter-dependent Dirichlet-multinomial log-likelihood.
- `multinomial_loglik`: multinomial log-likelihood used by the null test.
- `score_function`, `observed_fim`, `expected_fim`, `theta_fim`: numerical
  building blocks.
- `estimate_profile_loglik`, `grid_profile`, `adaptive_grid_profile`: fixed-
  theta profile-likelihood routines.
- `fit_equal_theta`: simultaneous fit of several count tables under one common
  theta.
- `random_dirichlet`: Dirichlet random generation.
- `sim_pop_sizes`, `sim_pop_equal_n`: Dirichlet-multinomial population
  simulation.
- `null_test`: simulation under H0: theta=0.
- `seed_rng`: deterministic seeding of the Fortran RNG used by this port.

See `API_COVERAGE.md` and `PORTING_NOTES.md` for the R-to-Fortran mapping and
intentional interface differences.

## Example

```fortran
program example
    use dirmult_fortran
    implicit none
    integer :: x(9,3)
    type(dirmult_fit_type) :: fit

    x = reshape([ &
        25,22,3,4,2,3,20,5,5, &
        3,6,24,23,4,5,5,20,5, &
        2,2,3,3,24,22,5,5,20], shape(x))

    call fit_dirmult(x, fit, epsilon=1.0e-10_dp, trace=.false.)
    print *, fit%pi
    print *, fit%theta
end program example
```

For this example the fitted theta is approximately `0.239353247`.
