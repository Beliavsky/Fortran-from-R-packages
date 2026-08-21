# Validation

The package is validated with GNU Fortran using strict runtime checking:

```text
-std=f2008
-O0 -g
-Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface
-fcheck=all
-ffpe-trap=invalid,zero,overflow
```

The v0.2.0 suite has 13 test programs covering:

- constant and stepped population rate tables
- expected population survival and transformed times
- `expprep2` consistency
- Pohar-Perme/Ederer/Hakulinen estimators
- additive Aalen increments
- direct additive relative-survival maximum likelihood
- `rsadd` EM, grouped binomial, and grouped Poisson paths
- `residuals.rsadd`-style residuals, including tied events
- `rs.br` and `rs.zph` diagnostics
- exact-boundary and predictable `epa` kernels
- HLD/HMD rate-table parsers
- years-difference, YL2013/YL2017 variance/bootstrap formulas
- follow-up splitting and rate-table utilities
- crude mortality (`cmp.rel`)
- grouped `rs.diff`
- analytical `nessie` expected counts

A clean direct build also runs `example/demo_relsurv.f90`.

Exact floating-point equality is intentionally retained in a few event-time/tie
comparisons because it mirrors the upstream C/R event-grid semantics. These
comparisons generate `-Wcompare-reals` warnings but do not fail the strict build.
