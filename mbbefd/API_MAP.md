# API map

## MBBEFD

| R | Fortran |
|---|---|
| `dmbbefd` | `dmbbefd` |
| `pmbbefd` | `pmbbefd` |
| `qmbbefd` | `qmbbefd` |
| `rmbbefd` | `rmbbefd` |
| `ecmbbefd` | `ecmbbefd` |
| `mmbbefd` | `mmbbefd` |
| `tlmbbefd` | `tlmbbefd` |
| `dMBBEFD` | `dmbbefd_gb` |
| `pMBBEFD` | `pmbbefd_gb` |
| `qMBBEFD` | `qmbbefd_gb` |
| `rMBBEFD` | `rmbbefd_gb` |
| `ecMBBEFD` | `ecmbbefd_gb` |
| `mMBBEFD` | `mmbbefd_gb` |
| `tlMBBEFD` | `tlmbbefd_gb` |
| `g2a` | `g2a` |
| `swissRe` | `swiss_re` |

The `_gb` suffix is required because Fortran is case-insensitive.

## Fitting helper densities

| R | Fortran |
|---|---|
| `dmbbefd1` | `dmbbefd1` |
| `dmbbefd2` | `dmbbefd2` |
| `dMBBEFD1` | `dmbbefd_gb1` |
| `dMBBEFD2` | `dmbbefd_gb2` |
| `dgbeta1` | `dgbeta1` |

## Other distributions

The following keep the same lower-case names:

- `dstpareto`, `pstpareto`, `qstpareto`, `rstpareto`, `ecstpareto`, `mstpareto`
- `dgbeta`, `pgbeta`, `qgbeta`, `rgbeta`, `ecgbeta`, `mgbeta`
- `ecunif`, `ecbeta`
- `doiunif`, `poiunif`, `qoiunif`, `roiunif`, `ecoiunif`, `moiunif`, `tloiunif`
- `doibeta`, `poibeta`, `qoibeta`, `roibeta`, `ecoibeta`, `moibeta`, `tloibeta`
- `doistpareto`, `poistpareto`, `qoistpareto`, `roistpareto`, `ecoistpareto`, `moistpareto`, `tloistpareto`
- `doigbeta`, `poigbeta`, `qoigbeta`, `roigbeta`, `ecoigbeta`, `moigbeta`, `tloigbeta`

The generic R helpers `doifun`, `poifun`, `qoifun`, `roifun`, `ecoifun`,
`moifun`, and `tloifun` are also provided. In Fortran they accept procedure
callbacks plus a real parameter vector instead of R's `...` argument.

## Empirical and fitting API

| R | Fortran |
|---|---|
| `etl` | `etl` |
| `eecf` | `make_eecf` returning `type(eecf_t)` |
| empirical EECF closure call | `eecf_t%evaluate(d)` |
| `fitDR` | `fit_dr` |
| `bootDR` | `boot_dr` |

`fit_dr` accepts distribution strings `oiunif`, `oistpareto`, `oibeta`,
`oigbeta`, `mbbefd`, and either `MBBEFD` or canonical `mbbefd_gb` for the
second parameterization.

## Intentionally omitted

- `eccomp`: graphics/comparison plotting
- S3 `print`, `summary`, `plot`, and `lines` methods for `eecf`
- Rcpp registration/glue and debug-only C++ entry points
