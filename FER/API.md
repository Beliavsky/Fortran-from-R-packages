# API

Use the umbrella module:

```fortran
use fer
```

The double-precision kind is `dp`. Option direction is represented by `cp=1`
for calls and `cp=-1` for puts.

## Vanilla options

- `bachelier_price(strike, forward, texp, sigma, cp [, df])`
- `bachelier_impvol(price, strike, forward, texp, cp [, df])`
- `black_scholes_price(strike, forward, texp, sigma, cp [, df])`
- `black_scholes_impvol(price, strike, forward, texp, cp [, df])`

## CEV

- `cev_price(strike, forward, texp, sigma, beta, cp [, df])`
- `cev_mass_zero(forward, texp, sigma, beta)`

## SABR and NSVh

- `sabr_hagan_2002(strike, forward, texp, sigma, vov, rho, beta)`
- `sabr_hagan_price(strike, forward, texp, sigma, vov, rho, beta, cp, df)`
- `nsvh1_choi_2019(strike, forward, texp, sigma, vov, rho, cp, df)`

## Exchange and spread options

- `switch_margrabe(forward1, forward2, texp, sigma1, sigma2, corr, cp, df)`
- `spread_kirk(strike, forward1, forward2, texp, sigma1, sigma2, corr, cp, df)`
- `spread_bjerksund_2014(strike, forward1, forward2, texp, sigma1, sigma2, corr, cp, df)`
- `spread_bachelier(strike, forward1, forward2, texp, sigma1, sigma2, corr, cp, df)`
