# blmodel-fortran

A dependency-free modern Fortran translation of the numerical core of the R
package `BLModel` 1.0.2.

The original package computes Black-Litterman posterior scenario distributions
from a discrete prior and a continuous investor-view distribution. This port
uses typed arrays, derived result types, and Fortran procedure callbacks while
preserving the original GPL version 3 license and attribution.

## Features

- Weighted discrete means and covariance matrices
- Inverse mean-risk equilibrium returns for CVaR, deviation CVaR, lower
  semi-absolute deviation, and mean absolute deviation
- Elliptical/variance equilibrium returns
- Normal, multivariate Student-t, and multivariate power-exponential view
  densities
- Diagonal or full view covariance matrices
- Posterior scenario reweighting
- Complete high-level Black-Litterman workflow
- Custom continuous view distributions through a procedure callback
- No external BLAS, LAPACK, R, or statistics-library dependency

## Build with FPM

```text
fpm build
fpm test
fpm run blmodel_demo
fpm run --example basic_black_litterman
fpm run --example view_distributions
```

A direct GNU Fortran validation script is also included:

```text
scripts/validate.sh
```

On Windows:

```bat
scripts\validate.bat
```

## Main interface

```fortran
use blmodel, only : dp, posterior_result, observ_normal, bl_post_distribution

type(posterior_result) :: result

result = bl_post_distribution(returns, probabilities, returns_freq=12.0_dp, &
  prior_type='elliptic', market_portfolio=weights, sharpe_ratio=0.5_dp, &
  pick=pick, q=annual_views, tau=0.05_dp, risk='MAD', alpha=0.95_dp, &
  view_density=observ_normal, view_covariance_type='diag')

if (.not. result%ok) error stop result%message
```

The result contains:

- `returns`: prior scenarios recentered at the equilibrium mean
- `probabilities`: normalized posterior scenario probabilities
- `equilibrium_returns`: equilibrium per-period expected returns
- `view_covariance`: covariance passed to the view-density callback
- `ok` and `message`: validation status

## Prior types

The R interface used `prior_type = "elliptic"` or `NULL`. Fortran uses explicit
strings:

- `elliptic` or `elliptical`
- `general`, `none`, `nonelliptic`, or `non-elliptic`

For a general prior, the risk argument may be `CVAR`, `DCVAR`, `LSAD`, or `MAD`.
As in the upstream code, CVAR and DCVAR share the same inverse-optimization
calculation, and LSAD and MAD share the same calculation.

## Custom view distributions

A callback must implement the public `view_density_interface`:

```fortran
subroutine my_density(points, q, covmat, params, density, info)
  use blmodel, only : dp
  real(dp), intent(in) :: points(:,:), q(:), covmat(:,:)
  real(dp), intent(in), optional :: params(:)
  real(dp), allocatable, intent(out) :: density(:)
  integer, intent(out) :: info
end subroutine my_density
```

`points` has shape `(number_of_views, number_of_scenarios)`, matching the actual
orientation used by the upstream R implementation. Each column is one point.
The optional `view_params` vector can carry distribution-specific parameters.
The built-in Student-t routine uses `view_params(1)` as degrees of freedom; the
power-exponential routine uses it as beta.

## Numerical behavior

The port normalizes input probabilities before computing moments and posterior
weights. The upstream package assumes they already sum to one. This produces
the same result for valid upstream inputs while making scale-equivalent weights
safe.

Distribution densities use Cholesky factors and log normalization constants
rather than explicit inverses and determinants. Covariance matrices must be
symmetric positive definite.

## Project structure

- `src/`: translated Fortran library
- `test/`: deterministic numerical and validation tests
- `app/`: demonstration program
- `example/`: focused examples
- `original/`: unmodified upstream package
- `provenance/`: supplied archive and checksum manifests
- `COVERAGE.md`: routine-by-routine mapping
- `PORTING_NOTES.md`: translation decisions
- `VALIDATION.md`: compiler and test record

## License

GPL-3.0-only. See `LICENSE` and `NOTICE`.
