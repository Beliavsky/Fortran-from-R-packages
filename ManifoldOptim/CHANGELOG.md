# Changelog

## 0.2.0

- Added constructed Stiefel ParamSet=2 retraction.
- Added Sphere parameter-set retractions, parallel translation, differentiated
  retraction, and locking Beta behavior.
- Added intrinsic-coordinate Stiefel/Grassmann transports.
- Added LowRank quotient projection, metric scaling, and D-aware transport.
- Replaced dense LRTRSR1 with compact limited-memory SR1.
- Added Armijo, weak Wolfe, strong Wolfe, exact, and custom line searches.
- Added RCG beta variants and changed the default to Hestenes-Stiefel.
- Added manifold Beta scaling and ROPTLIB-style quasi-Newton safeguards.
- Added full Broyden-family rank-one term and a configurable Phi extension.
- Added explicit cotangent actions for QF Stiefel/Grassmann and Sphere.
- Added Stiefel/Grassmann Euclidean-Hessian to Riemannian-Hessian corrections.
- Added architecture parity regression tests.

## 0.1.0

Initial standalone modern Fortran/FPM translation.
