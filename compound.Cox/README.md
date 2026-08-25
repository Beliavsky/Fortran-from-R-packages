# compound.Cox-fortran

**Official CRAN title:** Univariate Feature Selection and Compound Covariate for Predicting Survival, Including Copula-Based Analyses for Dependent Censoring

Modern Fortran/FPM translation of the numerical core of the R package
`compound.Cox` 3.33 by Takeshi Emura, Hsuan-Yu Chen, Shigeyuki Matsui and
Yi-Hau Chen.

The upstream package performs univariate feature screening and compound-covariate
survival prediction under ordinary Cox regression and Clayton-copula dependent
censoring. It also contains copula-graphic survival estimators and factorial
survival inference.

## Implemented numerical API

The public umbrella module is:

```fortran
use compound_cox
```

Main procedures:

- `cg_clayton`, `cg_frank`, `cg_gumbel`
- `cg_test`
- `uni_score`, `uni_wald`, `uni_selection`
- `compound_reg`
- `depend_cox_reg`, `cindex_cv`, `depend_cox_reg_cv`
- `surv_factorial`
- `x_pathway`, `x_tag`

Result structures are represented by Fortran derived types such as
`cg_result`, `compound_result`, `depend_cox_result`, `selection_result`, and
`factorial_result`.

## Dependencies

The supplied dependency translations are included under `vendor/`:

- `numDeriv-fortran` -- used for final likelihood Hessians in variance estimation.
- `survival-fortran-core` -- a compact FPM package assembled from the supplied
  survival translation, providing Cox fitting and concordance calculations.

The R package's `MASS::ginv` dependency is replaced by an internal symmetric
Jacobi-eigendecomposition pseudoinverse, so no MASS-equivalent library is needed.

## Build

```text
fpm test
fpm run --example demo_compound_cox
```

The project-owned source has also been tested directly with GNU Fortran 2018,
`-Wall -Wextra -Werror -fcheck=all`.

## Scope

Plotting, R data-frame/formula/S3 infrastructure, and the packaged `Lung`/`PBC`
R datasets are not converted. The underlying numerical methods are implemented.
Cross-validation uses the same deterministic sequential folds as upstream by
default; callers can shuffle rows before invoking the Fortran routines when the
upstream `randomize=TRUE` behavior is desired.

See `API_MAP.md` and `PORTING_NOTES.md` for details.
