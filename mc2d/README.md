# mc2d-fortran

Modern Fortran translation of the computational core of the R package
**mc2d 0.2.2** (Tools for Two-Dimensional Monte-Carlo Simulations).

The port targets Fortran 2018 and the Fortran Package Manager (FPM). Plotting,
ggplot2/ggpubr graphics, and R print-method infrastructure are intentionally
excluded. Numerical and simulation functionality is retained in typed Fortran
APIs.

## Major translated areas

- Bernoulli, generalized/scaled beta, subjective beta, bounded-lognormal,
  PERT, triangular, minimum-quantile-information, continuous empirical, and
  discrete empirical distributions.
- Dirichlet, multinomial, and multivariate-normal densities and random
  generation.
- Truncated inverse-CDF sampling and Latin hypercube sampling.
- `mcnode` types `0`, `V`, `U`, and `VU`, including broadcasting arithmetic,
  dimension control, variate extraction/combination, `pmin`/`pmax` analogues,
  and `outm` metadata.
- `mc` collections, summaries, per-variate multivariate summaries/quantiles,
  uncertainty/variability ratios, `mcapply`-style reductions, and running
  convergence statistics.
- Stochastic-node generation through Fortran procedure callbacks.
- Probability-tree sampling.
- Iman-Conover rank-correlation induction (`cornode`) for matrices and
  compatible `mcnode` arrays.
- Variability and uncertainty tornado correlation calculations, including
  multivariate `outm="each"` VU nodes.
- `mcmodel` / `mcmodelcut` evaluation through strongly typed procedure
  callbacks in place of R expression capture and dynamic evaluation.
  `evalmccut_reduce` provides the memory-saving per-uncertainty-column loop with
  a typed fixed-size reducer.

The supplied `mvtnorm-fortran` package is vendored as an FPM path dependency
and is used for multivariate-normal operations and linear algebra.


## Deliberate API adaptations

R vector recycling and S3/dynamic expression dispatch are expressed as typed
Fortran arrays, generic interfaces, derived types, and procedure callbacks.
The upstream ability to place a vector of arbitrary R function names in `outm`
is not reproduced dynamically: use `outm="each"`/`"none"`, or explicitly reduce
multivariate values in Fortran. Graphics and print-only methods are omitted as
requested. See `API_COVERAGE.md` and `PORTING_NOTES.md` for details.

## Build

```text
fpm build
fpm test
fpm run --example basic_mc2d
```

FPM was not installed in the translation environment, so validation was
performed by compiling the same source and test programs directly with GNU
Fortran 14.2 using Fortran 2018 mode and runtime bounds/checking enabled.
No nonstandard free-form line-length option is required.

## Quick API example

```fortran
use mc2d, only : dp, mcnode, mcdata, operator(+)
implicit none

type(mcnode) :: variability, uncertainty, total
integer :: i
real(dp) :: v(100), u(50)

do i=1,100
  v(i)=real(i,dp)/100.0_dp
end do
do i=1,50
  u(i)=0.1_dp*real(i,dp)/50.0_dp
end do

variability=mcdata(v,type='V',nsv=100)
uncertainty=mcdata(u,type='U',nsu=50)
total=variability+uncertainty

print *, shape(total%value) ! 100 50 1, type VU
```

## Licensing and provenance

The upstream `mc2d` package declares `GPL (>= 2)`, so this translation is
licensed under **GPL-2.0-or-later**. See `LICENSE`, `NOTICE`, and `provenance/`.
The vendored `mvtnorm-fortran` subtree retains its supplied license and notices.
