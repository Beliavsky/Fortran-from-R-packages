# R to Fortran API map

| Exported R function | Converted Fortran routine or API |
|---|---|
| `autcor` | `timsac::autocorrelation`; raw `autcorf` |
| `fftcor` | `fftcorf` |
| `mulcor` | `timsac::multivariate_correlation`; raw `mulcorf` |
| `auspec` | `timsac::power_spectrum`; raw `autcorf` plus `auspecf` |
| `mulspe` | `mulspef` |
| `sglfre` | `sglfref` |
| `mulfrf` | `mulfrff` |
| `fpeaut` | `fpeautf` |
| `fpec` | `fpec7` |
| `mulnos` | `mulnosf` |
| `raspec` | `raspecf` |
| `mulrsp` | `mulrspf` |
| `optdes` | `optdesf` |
| `optsim` | `optsimf` |
| `wnoise` | `timsac::white_noise`; raw `wnoisef` |
| `autoarmafit` | `autarm` |
| `armafit` | `autarm` |
| `bispec` | `bispecf` |
| `canarm` | `canarmf` |
| `canoca` | `canocaf` |
| `covgen` | `covgenf` |
| `markov` | `markovf` |
| `nonst` | `nonstf` |
| `prdctr` | `prdctrf` |
| `simcon` | `simconf` |
| `thirmo` | `thirmof` |
| `mfilter` | `timsac::matrix_filter` |
| `blocar` | `blocarf` |
| `blomar` | `blomarf` |
| `bsubst` | `bsubstf` |
| `exsar` | `exsarf` |
| `mlocar` | `mlocarf` |
| `mlomar` | `mlomarf` |
| `mulbar` | `mulbarf` |
| `mulmar` | `mulmarf` |
| `perars` | `perarsf` |
| `unibar` | `unibarf` |
| `unimar` | `unimarf` |
| `xsarma` | `xsarmaf` |
| `baysea` | `bayseaf` |
| `decomp` | `decompf` and `spgrh` |

R plotting methods, S3 print methods, time-series attributes, and graphics are
not numerical algorithms and are not reproduced. Numerical arrays returned by
the raw routines correspond to the arrays allocated in the original R
wrappers.
