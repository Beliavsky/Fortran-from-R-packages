# Porting notes

## Data layout

All public matrix APIs use observations in rows and variables or predictors in columns. For regression, add an explicit all-ones column when an intercept is required. Logistic compatibility procedures add an intercept by default, matching the R functions.

## Numerical dependencies

LAPACK and BLAS are required by the vendored robustbase implementation. The rrcov dependency supplies robust covariance/PCA primitives, probability distribution functions, random directions, and auxiliary linear algebra.

## Reproducibility

Randomized robust covariance starts accept integer seeds. Defaults are fixed, so repeated calls are reproducible unless a different seed is supplied.

## Status and convergence

Procedures return result types rather than raising R conditions. A result can contain useful estimates with `status = robstattm_no_convergence`; callers should inspect both `status` and `converged`. Invalid dimensions and singular systems receive distinct status codes.

## Differences from R object semantics

The Fortran result types hold only numerical outputs. They do not retain formulas, terms, row names, column names, factor contrasts, calls, model frames, or R environments. Column-level RFPE selection therefore cannot enforce R formula hierarchy rules automatically.

## License combination

RobStatTM is GPL-3.0-or-later. The robustbase translation is GPL-2.0-or-later, which permits use under GPL version 3. The rrcov translation is GPL-3.0-or-later. The distributed combined work is GPL-3.0-or-later.
