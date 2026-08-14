# anMC-fortran

Modern Fortran/FPM translation of the computational code in the R package
`anMC` 0.2.5 by Dario Azzimonti.

The library estimates high-dimensional Gaussian orthant probabilities using
the low-dimensional active-set decomposition and Monte Carlo / asymmetric
nested Monte Carlo (ANMC) bias correction described by Azzimonti and
Ginsbourger (2018).  It also implements conservative excursion-set estimates.

## Implemented API

The umbrella module is `anmc`.

- `proba_max` - estimate `P(max(X) > threshold)`.
- `proba_min` - estimate `P(min(X) < threshold)`.
- `anmc_gauss` - asymmetric nested Monte Carlo remainder estimator.
- `mc_gauss` - ordinary Monte Carlo remainder estimator.
- `conservative_estimate` - conservative Gaussian excursion-set estimate.
- `select_active_dims` - active-dimension heuristics 0 through 5.
- `select_q_dims` - sequential selection of the active-set size.
- `mvrnorm_arma` - multivariate-normal simulation, including Cholesky-input mode.
- `trmvrnorm_rej_cpp` - rejection sampler for a truncated multivariate normal.
- `chronotime_ns` - high-resolution elapsed-time counter corresponding to
  upstream `get_chronotime`.

Fortran-specific derived types (`probability_estimate`, `mc_result`,
`conservative_result`, `simulation_control`, and others) replace R lists and
attributes.

## Gaussian rectangle probabilities

The R package depends on `mvtnorm`.  A previously supplied Fortran translation
of `mvtnorm` is GPL-2.0-only, whereas `anMC` is GPL-3.  This release therefore
does not statically combine the two code bases.  Instead, the narrow
multivariate-normal rectangle-probability functionality needed by `anMC` is
implemented natively in `anmc_math` under GPL-3 using a Genz conditional
transformation with randomized Halton quasi-Monte Carlo.

The supplied `mvtnorm-fortran` package was retained as an independent
validation reference, not as a linked dependency.  The resulting FPM project
is self-contained and has no R, Rcpp, Armadillo, or external numerical-library
dependency.

## Example

```fortran
use anmc

integer, parameter :: n = 12
real(dp) :: mu(n), sigma(n,n), design(n,1)
type(probability_estimate) :: ans
integer :: i

mu = 0.0_dp
sigma = 0.5_dp
do i = 1, n
  sigma(i,i) = 1.0_dp
  design(i,1) = real(i-1,dp) / real(n-1,dp)
end do

call seed_fortran_rng(1234)
ans = proba_max(0.05_dp, 0.0_dp, mu, sigma, design, q=4, &
                method=0, algo='ANMC')
print *, ans%probability, ans%variance
```

See `example/equicorrelated_orthant.f90` for a complete program.

## Build

With FPM:

```text
fpm build
fpm test
fpm run --example equicorrelated_orthant
```

The release was also compiled directly with GNU Fortran 14.2 using:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

## Porting choices

R S3/list plumbing, package namespace code, garbage-collection calls, console
formatting, and plotting statements that occur only in examples are not
translated.  Numerical algorithms, timing logic, active-set selection,
Gaussian simulation, probability decomposition, bias correction, and
conservative-set computation are translated.

`simulation_control` supplies explicit safety limits for outer/inner sample
counts and rejection batches.  These avoid machine-dependent accidental
multi-gigabyte allocations when the Fortran implementation is substantially
faster than the original R/C++ timing calibration.  Defaults are intentionally
large enough for normal use and can be raised by the caller.

See `API_MAPPING.md`, `ALGORITHM_NOTES.md`, and `VALIDATION.md` for details.

## License

The translated `anMC` code is distributed under GPL-3.0-only, matching the
license declared by the supplied upstream package.  The complete GPL version 3
text is in `LICENSE`.  Upstream metadata is retained in `upstream/`.
