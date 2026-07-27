# Porting notes

## Source and license

The supplied repository identifies itself as `MarkowitzR` version
1.0.2.0002. `DESCRIPTION` declares `LGPL-3`; every computational R source file
states LGPL version 3 or, at the user's option, any later version. The Fortran
port therefore uses SPDX identifier `LGPL-3.0-or-later`.

## Data layout

R stores matrices in column-major order, as does Fortran. The port preserves
`matrixcalc::vech` ordering: each lower-triangular column is emitted from its
diagonal entry downward. This matters for `mu`, `covariance`, `w_indices`, and
`w_covariance` compatibility.

## R covariance dispatch

Upstream accepts an R function such as `vcov`, `vcovHC`, or `vcovHAC`. A
Fortran procedure cannot consume an R `lm` object, so the port provides:

- empirical covariance of the moment mean
- analytic covariance under multivariate normal returns
- Bartlett/Newey-West HAC covariance
- `moment_covariance_callback` for user-defined estimators

The callback receives the same underlying moment observations from which the R
intercept-only multivariate regression is formed.

## Analytic delta method

The inverse and projected-inverse covariance calculations use the analytic
Jacobian

```text
-L (P kron P) D
```

where `L` selects lower-triangular elements, `D` is the duplication matrix, and
`P` is the relevant inverse or projected precision. With hedging, the Jacobian
uses the difference of the two projected Kronecker products, matching upstream.

No finite-difference derivative is used.

## Weighting behavior

The upstream documentation says both returns and features are multiplied by
`weights`, but the R implementation multiplies only `feat`. When no feature
matrix is present, upstream silently ignores weights.

The port does not silently change that behavior. Its default is
`weights_upstream`. The documented transformation can be requested with
`weights_all_columns`. Tests verify that the two modes remain distinct.

The fitted intercept is not multiplied in either mode, matching the structure
of the upstream augmented moment.

## Constraint validation

The R source contains a FIXME to check that the row space of `Gmat` is spanned
by `Jmat`. The Fortran port performs this validation and returns a status error
when the condition fails. This prevents an invalid projected-difference model
from being reported as a valid hedge.

## Numerical safeguards

The port adds explicit checks for:

- inconsistent observation counts
- invalid weight lengths
- invalid constraint dimensions
- missing intercept and features simultaneously
- singular second moments
- singular constrained or hedging moments
- non-finite observations, removed by complete-case filtering

Covariance outputs are explicitly averaged across opposite triangles and then
mirrored. This makes symmetry exact across compilers and optimization levels.

## Dependencies

The R package depends on `matrixcalc` and `gtools` and can optionally use
`sandwich`. Their required computational pieces are implemented locally. The
Fortran library has no external dependencies.
