# chyper-fortran

Modern free-format Fortran translation of the computational code in the R
package `chyper` 0.3.1 by William Nickols.

The package implements the conditional hypergeometric distribution describing
the number of objects shared by all samples when each sample is drawn without
replacement from a population consisting of one common region plus a
sample-specific region.

## Main API

- `dchyper`, `dchyper_vec`
- `pchyper`, `pchyper_vec`
- `qchyper`, `qchyper_vec`
- `rchyper`, `rchyper_one`, `chyper_seed`
- `pvalchyper`
- `mle_s`, `mle_n`, `mle_m`
- `loglik_chyper`
- `chyper_probabilities`

All compiled code is pure Fortran 2018 in free-format `.f90` files. There is no
C/C++ runtime dependency.

Build with:

```text
fpm build
fpm test
```
