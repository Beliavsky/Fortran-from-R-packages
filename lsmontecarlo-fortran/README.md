# LSMonteCarlo modern Fortran port

This is a self-contained modern Fortran/FPM translation of the computational
parts of R package `LSMonteCarlo` 1.0 by Mikhail A. Beketov.

The library prices European and early-exercise put options using
Black-Scholes and Longstaff-Schwartz least-squares Monte Carlo methods. It
supports plain American puts, arithmetic-average Asian American puts, quanto
American puts, antithetic variates, a European-put control variate, geometric
Brownian-motion simulation, and volatility/strike price surfaces.

## Build and test

```text
fpm build
fpm test
fpm run lsmontecarlo_demo
fpm run --example basic_options
fpm run --example quanto_surface
```

The project has no external library dependencies. It requires a Fortran 2018
compiler supported by FPM.

A direct GNU Fortran validation script is also provided:

```text
bash scripts/validate_gfortran.sh
```

On Windows with GNU Fortran available in `PATH`:

```text
scripts\validate_gfortran.bat
```

## Main API

```fortran
use lsmontecarlo, only : dp, option_result
use lsmontecarlo, only : amer_put_lsm_av, eu_put_bs

type(option_result) :: result
real(dp) :: european

european = eu_put_bs(100.0_dp, 0.20_dp, 105.0_dp, 0.05_dp, 0.0_dp, 1.0_dp)
result = amer_put_lsm_av(100.0_dp, 0.20_dp, 10000, 50, 105.0_dp, &
    0.05_dp, 0.0_dp, 1.0_dp, 12345)

print *, european
print *, result%price, result%standard_error
```

Public pricing routines have both concise compatibility-style names and more
descriptive generic aliases:

- `amer_put_lsm` / `american_put_lsmc`
- `amer_put_lsm_av` / `american_put_lsmc_antithetic`
- `amer_put_lsm_cv` / `american_put_lsmc_control`
- `asian_amer_put_lsm` / `asian_american_put_lsmc`
- `quanto_amer_put_lsm` / `quanto_american_put_lsmc`
- `quanto_amer_put_lsm_av` / `quanto_american_put_lsmc_antithetic`

The original argument order and defaults are retained where practical. An
optional final integer `seed` argument was added to simulation-based routines.
Keyword arguments are recommended.

## Numerical implementation

The early-exercise routines use standard backward Longstaff-Schwartz
induction. Continuation values are estimated only on in-the-money paths.
Polynomial regressions use scaled, column-pivoted QR and automatically reduce
the basis rank when columns are dependent. This matters for degenerate cases,
such as a deterministic quanto multiplier.

Antithetic estimators report standard errors from paired path averages. The
control-variate estimator applies the European-put adjustment path by path.

## License and provenance

The original package declares `GPL-3`; this port therefore uses
`GPL-3.0-only`. The complete GPL version 3 text is in `LICENSE`. Every Fortran
source carries an SPDX identifier and original-package attribution.

The supplied package is retained unchanged under `original/LSMonteCarlo-1.0`.
Checksum manifests are under `provenance`.

See also:

- `COVERAGE.md` for the R-to-Fortran routine map
- `PORTING_NOTES.md` for behavioral and numerical details
- `VALIDATION.md` for the validation record
