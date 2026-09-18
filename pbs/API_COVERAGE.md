# API coverage

Coverage basis: distinct exported computational R functions and explicitly registered computational S3 methods in upstream `NAMESPACE`.

Package status: **substantial**.

Mapped: **2 of 2 (100%)**.

| R function | R source | Fortran module | Public Fortran procedure(s) | Status | Notes |
|---|---|---|---|---|---|
| `pbs` | `upstream/R/pbs.R` | `pbs_api` | `pbs` | substantial | Periodic and ordinary B-spline basis generation, type-7 quantile knot selection, missing rows, intercept handling, and nonperiodic Taylor extrapolation are translated. R attributes/classes/warnings are represented by a typed result and status/message interface. |
| `predict.pbs` | `upstream/R/pbs.R` | `pbs_api` | `predict_pbs` | substantial | Reuses the stored spline specification for new predictor values. R S3 dispatch and `...` are omitted. |

## Excluded upstream interface helpers

- `makepredictcall.pbs` is registered as an S3 method but only modifies an R language call used by formula/model-frame machinery. It is presentation/interface infrastructure rather than a numerical computation and is excluded from the coverage denominator under the repository coverage conventions.

## Material compatibility notes

- R matrix attributes, row/column names, `"pbs"`/`"basis"` classes, and S3 dispatch are not reproduced. Equivalent numerical metadata is stored in `type(pbs_basis)`.
- R warning emission is not reproduced verbatim. Invalid periodic boundary usage is returned through `status` and `message`; ordinary extrapolation is computed numerically as upstream does.
- User-supplied internal knots are retained in the returned metadata in their supplied order, while a sorted copy is used internally for basis construction.
