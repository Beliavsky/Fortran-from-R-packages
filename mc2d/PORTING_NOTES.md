# Porting notes

## Source baseline

The input package is mc2d 0.2.2, dated 2026-05-12. Its DESCRIPTION declares
`License: GPL (>= 2)`. The exact input archive and the original R sources used
for the translation are retained under `provenance/`.

## Design choices

### mcnode representation

An R `mcnode` is an array with dimensions `(nsv, nsu, nvariates)` and a type
attribute (`0`, `V`, `U`, or `VU`). The Fortran port uses a derived type with an
allocatable rank-3 real array plus an integer type code and `outm` metadata.
Arithmetic operators explicitly broadcast dimensions of size 1, reproducing
the computational role of the R methods without relying on implicit array
recycling.

### R dynamic evaluation

`mcstoc`, `mcmodel`, and `mcmodelcut` in R can receive function names or captured
expressions and then evaluate them dynamically. Modern Fortran has no direct
analogue. The port uses explicit procedure interfaces and procedure pointers.
This gives compile-time type checking and allows the same numerical workflows,
but callers express models as Fortran callbacks.

For the memory-saving `mcmodelcut` workflow, `evalmccut_reduce` executes an
optional setup callback once, a column callback once per uncertainty replicate,
and a reducer callback that returns a fixed-size statistic vector. This is the
typed Fortran counterpart of the upstream three-block loop. Arbitrary nested R
list results are intentionally not emulated. `evalmccut` remains available for
a full-model callback when looped reduction is not required.

### mvtnorm

The supplied `mvtnorm-fortran` translation is vendored under `vendor/` and used
rather than reimplementing its multivariate-normal and Cholesky facilities.
One behavioral normalization is made in the mc2d wrapper: the supplied
`dmvnorm_one` defaults to log density, while R `mc2d::dmultinormal` defaults to
ordinary density. The wrapper therefore passes `log_density=.false.` unless the
caller requests logarithms.

### Iman-Conover correlation

The upstream `cornode` constructs random normal-score permutations, ranks them,
multiplies those ranks by a Cholesky factor, ranks the result, and reorders each
input marginal. Because the initial normal scores are strictly ordered, their
ranks are exactly the sampled integer permutations. The Fortran implementation
uses those integer ranks directly, avoiding an unnecessary normal-quantile
calculation while preserving the upstream algorithm.


### Multivariate output metadata

The common upstream `outm="each"` and `outm="none"` cases are represented. R
also permits a vector of names of arbitrary R functions to be stored in `outm`
and looked up dynamically. The Fortran `mcnode` stores a single character
metadata value and does not perform dynamic function lookup. Multivariate
`tornado`/`tornadounc` support `outm="each"`; callers needing a custom reduction
should apply a Fortran reducer explicitly before analysis.

`node_summary_each` and `node_quantiles_each` provide one typed result per
variate for multivariate nodes.

### Numerical routines

The package uses an in-tree beta/noncentral-beta implementation built on the
regularized-beta routine from the supplied mvtnorm port. Random gamma, beta,
Poisson, binomial, and normal generation needed by mc2d are implemented in
`mc2d_random`.

## Validation

GNU Fortran 14.2 was used with `-std=f2018 -fcheck=all`. Tests cover:

- distribution identities and CDF/quantile inversion;
- PERT and triangular mean parameterizations;
- empirical distributions and duplicate aggregation;
- Dirichlet and multinomial constraints;
- multivariate-normal density and generation through mvtnorm-fortran;
- V/U/VU mcnode broadcasting and variate operations;
- LHS strata coverage and truncated sampling bounds;
- Iman-Conover rank correlation (target 0.8, observed about 0.8005 in the test);
- summaries (including per-variate multivariate summaries), tornado calculations,
  probability trees, full callback models, and looped `evalmccut_reduce` evaluation.
- multivariate VU tornado/tornadounc shape and correlation behavior.

FPM itself was not installed in the translation environment, so the FPM
manifest was not executed directly. The exact source/test units that FPM will
compile were compiled and linked manually without nonstandard line-length
flags.
