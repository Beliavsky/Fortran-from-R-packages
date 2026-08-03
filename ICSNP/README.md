# ICSNP-fortran

Modern Fortran translation of the computational core of the R package
**ICSNP 1.1-3**, *Tools for Multivariate Nonparametrics*.

The project is self-contained, uses Fortran Package Manager conventions, and
has no external numerical-library dependency.

## Implemented functionality

The public module `icsnp` provides:

- marginal rank, signed-rank, sign, and normal-score location tests;
- independent-component-model location and independence tests;
- one- and two-sample Hotelling T2 tests;
- Tyler-angle sign, rank, and normal-score location tests;
- spatial medians and spatial signs;
- Tyler, Duembgen, weighted Duembgen, symmetrized Huber, weighted Huber,
  HP1, and Hettmansperger-Randles estimators;
- Hodges-Lehmann and van der Waerden univariate location estimators;
- pairwise differences, sums, and products;
- spatial-rank and signed-rank kernels retained from the original C++ code.

R dots are replaced with underscores. For example, `tyler.shape` becomes
`tyler_shape`, and `HP.loc.test` becomes `HP_loc_test`.

## Build

```text
fpm build
fpm test
fpm run
fpm run --example location_example
fpm run --example shape_example
fpm run --example tests_example
```

GNU Fortran validation scripts are also included:

```text
./run_gfortran_tests.sh
run_gfortran_tests.bat
```

## Important porting differences

- R formula methods, S3 dispatch, `htest` objects, and `ics` object methods are
  replaced by explicit array-based interfaces and derived result types.
- The `ind_ictest` implementation uses a self-contained FOBI-style invariant
  coordinate transform. The R package delegates this step to the external
  `ICS` package.
- Simulation and permutation tests use a deterministic portable pseudo-random
  generator when a seed is supplied.
- Missing-value filtering is not implicit. Callers should provide finite
  arrays, making data-cleaning behavior explicit.

See `PORTING.md`, `API.md`, and `TRANSLATION_COVERAGE.md` for details.

## License

GPL-2.0-or-later, matching the original package metadata. Original source
material is retained under `original/`.
