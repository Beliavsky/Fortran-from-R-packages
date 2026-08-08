# Porting notes

## One-dimensional integration

The Gauss-Kronrod nodes and weights are generated directly from the supplied
`GaussKronrodNodesWeights.h`. The adaptive algorithm follows the package's
QUADPACK-style interval bisection and error rescaling.

The public infinite-interval transformations are preserved:

- `[a, +Inf)`: `x = a + (1-t)/t`;
- `(-Inf, b]`: `x = b - (1-t)/t`; and
- `(-Inf, +Inf)`: symmetric evaluation at `+x` and `-x`.

### Corrected 121-point rule defect

The supplied C++ integrator decides whether the embedded Gauss rule has a
center node from the parity of the enum value. This is incorrect for the
121-point Kronrod rule, whose embedded Gauss rule has 60 points and therefore
no center node. The Fortran implementation selects center-node rules from the
actual Gauss orders (7, 15, 25, 35, and 45).

Without this correction, smooth polynomial integrals have accurate values but
spurious error estimates near `1e-3` under rule 121.

## Multidimensional integration

The Fortran Cuhre module directly ports the rule construction, fully symmetric
point expansion, fourth-difference split-dimension selection, null-rule error
estimate, region bisection, and error correction used by the supplied Cuba
source. Since RcppNumerical calls Cuhre with `key=0`, the translated paths are:

- rule 13 in two dimensions;
- rule 11 in three dimensions; and
- rule 9 in four or more dimensions.

RcppNumerical sets Cuba's `LAST` flag, so the public result is the sum of the
latest region estimates. State files, parallel worker processes, multiple
integrand components, and verbose output are not part of the public wrapper
and are omitted.

### Corrected mixed-bound Jacobian defect

In the supplied `MFuncWithInfiniteBounds`, finite dimensions are linearly
mapped from `[0,1]` but their range factors are not included in the Jacobian
when another dimension has an infinite bound. The Fortran implementation
includes every finite-dimension range factor.

## Optimization

The root package uses two local FPM dependencies:

- `lbfgs`, a modern Fortran L-BFGS/libLBFGS translation; and
- `lbfgsb3`, a modern Fortran L-BFGS-B 3.0 translation.

The wrapper maps the RcppNumerical controls (`maxit`, `eps_f`, and `eps_g`) to
the closest native controls. The public callback contract is unchanged in
substance: one procedure returns both objective and gradient.

## R-specific code

Rcpp/Eigen maps, R warnings, external pointers, dynamic registration, and R
list construction are replaced by explicit arrays and derived result types.
The package contains no plotting code.
