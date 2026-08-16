# R to Fortran API map

All entries below are provided by `use sadists` unless marked omitted.

| R family | Fortran density | CDF | Quantile | RNG |
|---|---|---|---|---|
| doubly noncentral beta | `ddnbeta` | `pdnbeta` | `qdnbeta` | `rdnbeta` |
| doubly noncentral eta | `ddneta` | `pdneta` | `qdneta` | `rdneta` |
| doubly noncentral F | `ddnf` | `pdnf` | `qdnf` | `rdnf` |
| doubly noncentral t | `ddnt` | `pdnt` | `qdnt` | `rdnt` |
| K-prime | `dkprime` | `pkprime` | `qkprime` | `rkprime` |
| lambda-prime | `dlambdap` | `plambdap` | `qlambdap` | `rlambdap` |
| product chi-square powers | `dprodchisqpow` | `pprodchisqpow` | `qprodchisqpow` | `rprodchisqpow` |
| product doubly noncentral F | `dproddnf` | `pproddnf` | `qproddnf` | `rproddnf` |
| product normals | `dprodnormal` | `pprodnormal` | `qprodnormal` | `rprodnormal` |
| sum chi-square powers | `dsumchisqpow` | `psumchisqpow` | `qsumchisqpow` | `rsumchisqpow` |
| sum log chi-square | `dsumlogchisq` | `psumlogchisq` | `qsumlogchisq` | `rsumlogchisq` |
| upsilon | `dupsilon` | `pupsilon` | `qupsilon` | `rupsilon` |

`runExample` is omitted because it is a Shiny UI launcher.

Additional public computational helpers expose moment/cumulant machinery used
by the package, including `dnf_moments`, `dnt_moments`,
`chipow_cumulants`, `sumchisqpow_cumulants`, `sumlogchisq_cumulants`,
`lambdap_cumulants`, `kprime_cumulants`, `upsilon_cumulants`, and
`prodnormal_cumulants`.

## PDQutils integration

The v0.1 public helper names `edgeworth_pdf`, `edgeworth_cdf`,
`cornish_fisher_quantile`, `as269`, `moments_to_cumulants`, and
`cumulants_to_moments` remain available through `use sadists`. In v0.2.0 they
are compatibility views/wrappers over the vendored standalone
`pdqutils-fortran` implementation rather than independent algorithms.

