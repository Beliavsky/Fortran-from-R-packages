# Porting notes

## R object system

The original package exposes a large formula/S3 interface. Fortran has no
corresponding dynamic formula or method-dispatch runtime, so this port accepts
explicit response vectors, design matrices, grouping vectors, covariates, and
coordinates. Model results are derived types rather than S3 objects.

## GLS and LME likelihoods

`fit_gls` profiles the common residual scale and optimizes the selected ML or
REML likelihood. `fit_lme` constructs the marginal block covariance

```text
V_i = Z_i G Z_i^T + sigma^2 D_i R_i D_i
```

and optimizes its ML or REML likelihood. Fixed effects are generalized
least-squares estimates and random effects are Gaussian BLUPs.

The optimizer is a deterministic self-contained Nelder-Mead implementation.
This differs from the collection of R optimizers and native routines selectable
by the original package, so iteration paths and final low-order digits can
vary.

## Nonlinear models

`fit_gnls` uses a finite-difference Jacobian and a damped
Levenberg-Marquardt iteration. Supplied residual covariance structures are used,
but their parameters are not currently re-estimated during the nonlinear
iteration.

`fit_nlme` uses first-order iterative linearization. At each outer iteration it
builds fixed and random Jacobian design matrices and calls the translated linear
mixed-effects solver. This represents the main computational model, but it is
not a line-by-line reproduction of every PNLS and native optimization detail in
R `nlme`.

## Correlation structures

AR(1), CAR(1), compound symmetry, spatial kernels, and ARMA autocorrelations use
the same mathematical definitions as the original package. `COR_UNSTRUCTURED`
uses a normalized lower-Cholesky parameterization to guarantee positive
definiteness; this is numerically equivalent in scope to `corSymm`/`corNatural`
but not the same internal coordinate system.

## Random-effect structures

The original `pdSymm`, `pdLogChol`, and `pdNatural` classes are represented by a
stable log-Cholesky parameterization. `pdBlocked` and nested `reStruct` objects
are represented by assembling the desired block structure directly in `Z` or by
calling multiple models, rather than reproducing R container classes.

## Omitted infrastructure

Plotting, lattice panels, formula parsing, grouped-data classes, S3 printing and
replacement methods, R data-frame operations, package datasets, confidence-
interval presentation methods, and R-specific update/predict plumbing are not
compiled. Relevant original source files remain under `original/`.
