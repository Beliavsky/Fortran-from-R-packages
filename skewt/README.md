# skewt-fortran

Modern Fortran translation of the computational core of the R package
`skewt` (Robert King and Emily Anderson), implementing the Fernandez-Steel
skewed Student-t distribution.

## API

```fortran
use skewt, only : dp, dskt, pskt, qskt, rskt
```

- `dskt(x, df [, gamma])` density
- `pskt(x, df [, gamma])` CDF
- `qskt(p, df [, gamma])` quantile
- `rskt(x, df [, gamma])` random sample into array `x`

Scalar and rank-1 array forms are supplied for d/p/q.  `gamma=1` is the
ordinary Student-t distribution.  `df` and `gamma` must both be positive.

The implementation is self-contained.  Student-t probabilities use the
regularized incomplete beta function; quantiles are obtained by monotone
inversion; random generation follows the upstream inverse-transform design.

## Build

```text
fpm test
fpm run --example demo_skewt
```

No external libraries are required.

See `NOTICE.md` and `PORTING_NOTES.md` for provenance and translation notes.
