# API mapping

All 11 exported upstream numerical generators have Fortran counterparts:

- `generate_halton_faure_single` -> `generate_halton_faure_single`
- `generate_halton_random_single` -> `generate_halton_random_single`
- `generate_halton_faure_set` -> `generate_halton_faure_set`
- `generate_halton_random_set` -> `generate_halton_random_set`
- `generate_sobol_set` -> `generate_sobol_set`
- `generate_sobol_owen_set` -> `generate_sobol_owen_set`
- `generate_pj_set` -> `generate_pj_set`
- `generate_pmj_set` -> `generate_pmj_set`
- `generate_pmjbn_set` -> `generate_pmjbn_set`
- `generate_pmj02_set` -> `generate_pmj02_set`
- `generate_pmj02bn_set` -> `generate_pmj02bn_set`

The umbrella module also exposes low-level `sobol_single` and
`sobol_owen_single` functions.

For low-level/single-point routines, sequence index and dimension are zero
based, matching the C++ numerical kernels. Set routines use normal Fortran
matrix storage and return rows 1:n corresponding to sequence indices 0:n-1.

The supplied upstream C++ wrappers for the two Halton `*_single` routines call
`Halton_sampler::sample(i, dim)` even though the sampler signature is
`sample(dimension, index)`. The Fortran port deliberately corrects that
argument-order defect and follows the sampler's intended `(dimension,index)`
semantics.

Rcpp registration, R list/matrix conversion, roxygen documentation machinery,
and plotting examples are omitted.
