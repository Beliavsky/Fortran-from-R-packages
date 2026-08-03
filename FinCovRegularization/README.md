# FinCovRegularization for modern Fortran

A dependency-free modern Fortran translation of the computational routines in
**FinCovRegularization 1.1.0**, packaged for the Fortran Package Manager (FPM).
The original R package estimates and regularizes covariance matrices for
financial applications and constructs minimum-variance and risk-parity
portfolios.

## Features

- Squared Frobenius and spectral/operator matrix norms
- Banding, tapering, hard-thresholding, and soft-thresholding
- Deterministic repeated-split cross-validation for regularization parameters
- Minimum threshold search for positive definiteness
- Macroeconomic, fundamental OLS/WLS, and statistical factor covariance models
- Independence covariance operator matching the original R implementation
- Global minimum-variance portfolios with or without short sales
- Risk-parity optimization using the original package objective
- Explicit status codes and no required external numerical libraries

Plotting and R-specific S3 display methods are intentionally omitted. The
original R source, documentation, metadata, license, and data file are retained
under `original/`.

## Build with FPM

```text
fpm build
fpm test
fpm run
fpm run --example regularization_example
fpm run --example factor_models_example
fpm run --example portfolio_example
```

The package uses no FPM dependencies. It requires a Fortran 2018 compiler.

## Minimal example

```fortran
program covariance_example
   use fincovregularization
   implicit none

   real(dp) :: sigma(3,3), regularized(3,3), weights(3)
   integer :: status

   sigma = reshape([ &
      0.040_dp, 0.018_dp, 0.012_dp, &
      0.018_dp, 0.090_dp, 0.015_dp, &
      0.012_dp, 0.015_dp, 0.160_dp], [3,3])

   regularized = soft_thresholding(sigma, 0.01_dp, status)
   weights = gmvp(regularized, allow_short=.false., status=status)

   print '(3f12.6)', weights
end program covariance_example
```

## API naming

R dots are converted to underscores and names are lower case:

| R name | Fortran name |
|---|---|
| `F.norm2` | `f_norm2` |
| `O.norm2` | `o_norm2` |
| `FundamentalFactor.Cov` | `fundamental_factor_cov` |
| `MacroFactor.Cov` | `macro_factor_cov` |
| `StatFactor.Cov` | `stat_factor_cov` |
| `Ind.Cov` | `ind_cov` |
| `GMVP` | `gmvp` |
| `RiskParity` | `risk_parity` |
| `banding.cv` | `banding_cv` |
| `tapering.cv` | `tapering_cv` |
| `threshold.cv` | `threshold_cv` |
| `threshold.min` | `threshold_min` |

See [API.md](API.md) for argument and result details and [PORTING.md](PORTING.md)
for translation decisions.

## License

The original package declares GPL-2 and includes the GNU General Public License,
version 2. This translation is distributed under **GPL-2.0-only**. See
[LICENSE](LICENSE) and [NOTICE.md](NOTICE.md).
