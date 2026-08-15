# Delaporte-fortran

Standalone modern Fortran/FPM translation of the computational code from the R
package **Delaporte 8.4.3** by Avraham Adler.

The Delaporte distribution here uses the upstream parameterization and can be
viewed as a Poisson component with mean `lambda` plus a negative-binomial
component induced by a Gamma(`alpha`, scale=`beta`) mixed Poisson rate.
Consequently

- `E[X] = alpha*beta + lambda`
- `Var[X] = alpha*beta*(1 + beta) + lambda`

## Implemented computational API

- `ddelap` / `ddelap_vec`: probability mass function
- `pdelap` / `pdelap_vec`: cumulative distribution function
- `qdelap` / `qdelap_vec`: exact quantile by summation/inversion
- `rdelap` / `rdelap_vec`: random variates
  - `exact=.true.`: uniform variates followed by exact quantile inversion,
    matching the upstream compiled path
  - `exact=.false.`: Gamma-Poisson mixture generation, matching the upstream
    R fallback algorithm
- `qdelap_approx`: Monte-Carlo quantile approximation corresponding to
  upstream `qdelap(..., exact=FALSE)`; optional `nsim` lets callers override
  the upstream sample-size rule
- `momdelap`: method-of-moments estimates with skewness types 1, 2, and 3
- `seed_delaporte`: convenience wrapper around `random_seed`

Vector routines recycle parameter vectors in the same style as the R package.
Invalid non-positive parameters produce IEEE quiet NaNs. The scalar PMF returns
zero for non-integer observations, as upstream does.

The OpenMP thread getter/setter and R registration/C wrapper code are omitted:
they are interface/runtime plumbing rather than distribution algorithms. The
standalone library has no R dependency.

## Build with FPM

```text
fpm build
fpm test
fpm run --example basic
```

## Minimal use

```fortran
program demo
    use delaporte, only : dp, ddelap, pdelap, qdelap
    implicit none

    print *, ddelap(4.0_dp, 1.0_dp, 4.0_dp, 2.0_dp)
    print *, pdelap(4.0_dp, 1.0_dp, 4.0_dp, 2.0_dp)
    print *, qdelap(0.4_dp, 1.0_dp, 4.0_dp, 2.0_dp)
end program demo
```

## Licensing and provenance

BSD-2-Clause. See `LICENSE` and `NOTICE.md`. The `upstream/` directory contains
copies of the original computational Fortran sources plus DESCRIPTION/CITATION
for provenance. No plotting code exists in the upstream package.
