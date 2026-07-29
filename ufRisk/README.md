# ufRisk-fortran

A modern Fortran 2018/FPM translation of the computational algorithms in the R package **ufRisk 1.0.7**.

The library estimates one-step-ahead volatility, Value at Risk (VaR), and Expected Shortfall (ES) from six parametric or semiparametric volatility models and evaluates those forecasts with the original package's traffic-light, coverage, and loss-function procedures.

## Implemented models

- standard GARCH (`sGARCH`)
- exponential GARCH (`eGARCH`)
- asymmetric power ARCH (`apARCH`)
- fractionally integrated GARCH (`fiGARCH`)
- Log-GARCH through its ARMA representation (`lGARCH`)
- fractionally integrated Log-GARCH through its ARFIMA representation (`filGARCH`)

Normal and standardized Student-t innovations are supported. Each model may be run parametrically or with a local-polynomial nonparametric scale function.

## Implemented risk and validation algorithms

- rolling one-step conditional-volatility forecasts
- VaR at an independently selectable confidence level
- ES at an independently selectable confidence level
- VaR and ES traffic-light tests
- Kupiec unconditional-coverage test
- Christoffersen independence and conditional-coverage tests
- four regulatory, firm, Abad, and Feng loss functions

## Dependencies included in the source tree

This is a self-contained source distribution. It contains the numerical modules needed from the prior Fortran ports of:

- `rugarch` for sGARCH, eGARCH, APARCH, and FIGARCH
- `fracdiff` for ARFIMA estimation and FI-Log-GARCH
- `smoots` for local-polynomial scale estimation and ARMA fitting

The long-memory scale branch embeds the `tsmoothlm` workflow required by ufRisk: iterative local-polynomial smoothing, ARFIMA order selection by BIC, fractional-memory estimation, long-memory bandwidth inflation, and scale normalization. It is not a standalone translation of every `esemifar` function.

## Build

```sh
fpm build
fpm test
fpm run ufrisk_demo
```

A GNU Fortran fallback script is also provided:

```sh
./scripts/test_gfortran.sh
```

## Minimal use

```fortran
use ufrisk

type(ufrisk_options) :: options
type(ufrisk_result) :: result

options%model = ufrisk_model_sgarch
options%distribution = ufrisk_distribution_student
options%n_out = 250
result = varcast(prices, options)
```

See `API.md`, `PORTING.md`, and the programs under `app/`, `example/`, and `test/`.

## Scope

The R package's plotting, S3 printing, startup messages, bundled market datasets, and R list/data-frame presentation are omitted. All four exported computational entry points are translated.

## License

GPL-3.0-only. See `LICENSE` and `NOTICE.md`.
