# Translation coverage

## Included computational code

The port implements the numeric behavior represented by the exported rrcov
families:

- robust location/scatter estimation;
- robust and classical principal component analysis;
- robust and classical linear/quadratic discriminant analysis;
- Hotelling and Wilks multivariate tests;
- robust scale, outlyingness, matrix, ranking, and compositional helpers.

The roles of all five upstream native kernels are covered by free-form Fortran:

- `covOPW.c` -> `cov_ogk`;
- `fast-mve.c` -> `cov_mve`;
- `sest.c` -> `cov_sest`;
- `ds11.f` -> `cov_sde`;
- `fsada.f` -> grouped robust covariance through `cov_mwcd` and `linda_fit`.

## Deliberately omitted

- All plotting, biplots, score plots, distance plots, scree plots, ellipses,
  pair panels, and interactive point labeling.
- R formula/model-frame processing and missing-value class bookkeeping.
- S3/S4 class definitions, slots, accessors, `show`, `summary`, and R call
  reconstruction.
- R-specific random-number-state adapters, dynamic registration, and `.C`,
  `.Call`, or `.Fortran` interfaces.
- Data sets and serialized `.rda` loading APIs. The upstream data files remain
  in the provenance snapshot.
- Monte Carlo calibration closures returned by the R Wilks implementation.
  The Fortran API provides deterministic Bartlett and Rao approximations.

## Algorithmic fidelity

Classical covariance, PCA, LDA/QDA, Hotelling formulas, Wilks formulas, robust
univariate scales, OGK pairwise constructions, and concentration-step logic are
direct numerical translations or standard equivalent formulations.

Some upstream estimators delegate major work to `robustbase` or `pcaPP` rather
than containing the implementation in rrcov. To keep this package
self-contained, the Fortran versions implement documented algorithmic
counterparts. They are intended for robust numerical analysis and API coverage,
but bit-for-bit equality with every `robustbase`/`pcaPP` version is not claimed.
See `PORTING_NOTES.md`.
