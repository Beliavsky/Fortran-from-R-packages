# bivgeom-fortran

Modern, self-contained Fortran translation of the computational code in the
R package `bivgeom` 1.0 (Alessandro Barbiero), which implements Roy's
bivariate geometric distribution.

All compiled source is free-format Fortran 2018 (`.f90`). There is no C/C++
runtime code and no external numerical dependency.

## Computational coverage

- Roy bivariate geometric PMF and log-PMF
- joint CDF and joint survival function
- conditional CDF of `Y | X=x`
- conditional mean `E[Y | X=x]`
- local failure rates `lambda1` and `lambda2`
- Pearson correlation
- stress-strength reliability `P(X <= Y)`
- random generation
- empirical joint survival function
- negative log likelihood
- maximum-likelihood estimation
- least-squares estimation
- MMP and MM1--MM4 estimators

The module `bivgeom` exposes both explicit Fortran-style names such as
`dbivgeom_roy` and familiar compatibility names such as `dbivgeomroy`,
`fbivgeomroy`, `rbivgeomroy`, `estbivgeomroy`, and `minuslogroy`.

## Build

```text
fpm build
fpm test
fpm run --example roy_example
```

The release was validated directly with GNU Fortran using:

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all
```

## Example

```fortran
use, intrinsic :: iso_fortran_env, only : int64
use bivgeom
implicit none

integer, allocatable :: z(:, :)
type(bivgeom_fit) :: fit

call seed_rng(int(12345, int64))
call rbivgeom_roy(2000, 0.5_dp, 0.7_dp, 0.9_dp, z)
fit = fit_bivgeom_ml(z(:, 1), z(:, 2))
print *, fit%theta
```

See `TRANSLATION_NOTES.md` for the detailed R-to-Fortran mapping and numerical
differences.
