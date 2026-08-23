# PSDistr-fortran

Modern Fortran translation of the computational API of the R package **PSDistr** 0.0.1.

Implemented distribution families:

- Two-piece power normal: `dtppn`, `ptppn`, `qtppn`, `rtppn`
- Plasticizing component: `dpc`, `ppc`, `qpc`, `rpc`
- DS normal: `ddsn`, `pdsn`, `qdsn`, `rdsn`
- Expnormal: `den`, `pen`, `qen`, `ren`
- Sulewski plasticizing component: `dspc`, `pspc`, `qspc`, `rspc`
- Easily changeable kurtosis: `deck`, `peck`, `qeck`, `reck`

Scalar and rank-1 array interfaces are supplied for d/p/q functions. RNG routines fill a caller-provided array.

## Build

```sh
fpm test
fpm run --example demo_psdistr
```

The numerical core has no required external dependency. A copy of the supplied `pracma-fortran` translation is retained in `vendor/pracma-fortran` for provenance and comparison; PSDistr uses pracma upstream only for `nthroot`, for which this port uses an equivalent signed real-root helper.

## License

GPL-3.0-only, matching the upstream PSDistr DESCRIPTION. See `NOTICE.md` and `LICENSE`.
