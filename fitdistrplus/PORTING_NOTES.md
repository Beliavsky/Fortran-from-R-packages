# Porting notes

## Design

The R package discovers `d*`, `p*`, `q*`, `r*`, and moment functions by name.
Fortran uses `type(distribution_model)` with procedure-pointer components. This
keeps the algorithms generic while retaining compile-time interfaces.

Parameter domains are transformed to unconstrained coordinates before
Nelder-Mead optimization:

- finite lower and upper bounds: logistic transformation
- finite lower bound: shifted exponential
- finite upper bound: reflected exponential
- no bounds: identity

The MLE covariance matrix is the inverse numerical Hessian in transformed
coordinates, mapped back by the transformation Jacobian.

## Source-parity details

The MGE objective formulas preserve the scale used by `mgedist`, which differs
from the reporting scale used by `gofstat` for AD and CvM.

For weighted maximum spacing, the upstream implementation sorts observations
but does not explicitly reorder weights. This port reorders each weight with
its observation, which is the mathematically consistent interpretation.

The upstream uniform MLE closed form is retained.

## Deliberate omissions or replacements

- R plotting and graphical comparison functions
- S3 objects, formulas, data frames, calls, and expression capture
- R optimizer selection and arbitrary `custom.optim` functions
- R named-list `fix.arg` management
- parallel bootstrap execution
- asymptotic MME covariance and distribution-specific R test labels
- survival-package Kaplan-Meier and Turnbull object machinery

A custom Fortran distribution can represent a model with fixed parameters by
exposing only its free parameters through the callbacks.

## Random numbers

The library uses the compiler's `random_number` stream, Box-Muller normals,
Marsaglia-Tsang gamma generation, and standard mixture constructions. Seeded
runs are deterministic within a given compiler/runtime but are not bitwise
identical to R's RNG streams.
