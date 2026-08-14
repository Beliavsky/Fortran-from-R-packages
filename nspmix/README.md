# nspmix-fortran

Modern free-form Fortran/FPM translation of the computational core of R package
`nspmix` 2.0-0 (Yong Wang), for nonparametric and semiparametric mixture
maximum-likelihood estimation.

## Implemented numerical functionality

- finite discrete mixing distributions (`disc_dist`), sorting, normalization and collapse
- normal, Poisson, geometric and negative-binomial nonparametric mixture families
- common-variance random-effects problem (`CVPS`)
- grouped binomial logistic model with a nonparametric random intercept (`mlogit`)
- component log densities and analytic derivatives in support and structural parameters
- mixture log likelihood, mixture density and directional-gradient calculations
- hierarchical/constrained-Newton mass fitting (`hcnm`)
- support-search NPMLE (`cnm`)
- fixed-support proportion fitting (`cnm_proportions`)
- semiparametric compatibility entry points (`cnmms`, `cnmpl`, `cnmap`)
- normal/Poisson/geometric/negative-binomial mixture pdf/pmf and CDF routines
- random generation for the four basic univariate mixture families
- weighted histogram and CVP sufficient-statistic construction

The package is self-contained. The supplied `lsei-fortran-v0.1.0` sources are
vendored under `src/vendor_lsei` and used for equality-constrained,
nonnegative least-squares subproblems.

## Important implementation note

Upstream `nspmix` calls `lsei::pnnls(..., sum=delta)` in the CNM algorithms.
The supplied Fortran dependency's `pnnls_solve(..., sum_value=...)` returned a
zero vector on a direct simplex regression used during integration. To avoid
changing the dependency source, this port expresses the mathematically same
subproblem through `lsei_solve`: minimize `||A p-b||` subject to `sum(p)=delta`
and `p>=0`. This uses the supplied LSEI active-set constrained least-squares
solver and produced the expected likelihood optimum in randomized tests.

The upstream `cnmms`, `cnmpl` and `cnmap` routines use different BFGS/line-search
schedules. The v0.1.0 Fortran compatibility entry points share the same
profile/alternating semiparametric optimization engine while preserving the
same likelihood, constraints, support search and structural parameters. This
is the principal algorithmic implementation difference in v0.1.0.

## FPM

```toml
[fortran]
implicit-typing = false
implicit-external = false
source-form = "free"
```

No R runtime is required.

## Build

```text
fpm build
fpm test
fpm run --example poisson_mixture
```

## License

The upstream package is GPL (>= 2). The supplied LSEI translation is
GPL-2.0-or-later. This combined Fortran distribution is therefore distributed
under GPL-2.0-or-later. See `COPYING`, `UPSTREAM.md`, and the retained upstream
metadata.
