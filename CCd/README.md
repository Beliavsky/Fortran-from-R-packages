# CCd-fortran

Modern Fortran 2018/FPM translation of **CCd 1.1**, the Cauchy-Cacoullos
(discrete Cauchy) distribution package.

## Features

* `dcc`: PMF and log-PMF
* `pcc`: upstream-compatible CDF computation
* `qcc`: left quantiles
* `cc_mle`: free-location maximum likelihood fit
* `cc_mle0`: zero-location maximum likelihood fit
* `cc_reg`: Cauchy-Cacoullos regression with an intercept
* `loc0_test`: likelihood-ratio test for zero location
* standalone Nelder-Mead and golden-section optimizers
* no R or Rfast runtime dependency

## Example

```fortran
use ccd, only : dp, i8, dcc, pcc, qcc

print *, dcc(0.0_dp, 0.0_dp, 1.5_dp)
print *, pcc(2_i8, 0.0_dp, 1.5_dp)
print *, qcc(0.5_dp, 0.0_dp, 1.5_dp)
```

See `PORTING_NOTES.md` for an important upstream caveat concerning nonzero
location values in `pcc`/`qcc` and continuous-location likelihood fitting.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The package is GPL-2.0-or-later, matching upstream CCd.
