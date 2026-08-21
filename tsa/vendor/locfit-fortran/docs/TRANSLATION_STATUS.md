# Translation status: v0.1.0

The goal of v0.1.0 is a usable modern-Fortran local-fitting core rather than a
binary-compatible rewrite of the R/C implementation.  Computational behavior
was translated at the algorithm level and exposed through typed Fortran APIs.

## Substantially translated

- Kernel definitions, radial/product distances and kernel moments/convolutions.
- Standard polynomial/angular fitting basis.
- Family/link likelihood equations, including censored likelihood branches.
- Local weighted Newton likelihood solve and Gaussian weighted-least-squares
  specialization.
- Nearest-neighbor/fixed bandwidth neighborhoods.
- Automatic variable scaling.
- Coefficients, fitted means, local log likelihood, covariance and standard
  errors.
- Local derivatives.
- One-dimensional local density and ordinary KDE.
- KDE bandwidth-selection criteria from `band.c`.
- Residual and common model-selection criterion calculations.
- Selected high-level robust/quasi iterative procedures.
- Interpolation primitives from `ev_interp.c`.

## Preserved upstream but not yet given a native v0.1.0 API

The following specialized computational subsystems remain in the preserved
upstream tree and are candidates for later Fortran versions:

- adaptive local bandwidth selection from `lf_adap.c` and the full `regband`
  machinery;
- minimax smoothing weights from `minmax.c`;
- full local influence/hat-matrix and derivative-correction machinery from
  `lf_wdiag.c` and `lf_dercor.c`;
- parametric-component decomposition from `pcomp.c`;
- multidimensional density integration and specialized hazard integrators from
  `dens_int.c`, `dens_haz.c`, `dens_odi.c`, `m_isimp.c`, `m_isphr.c`, and
  `m_imont.c`;
- adaptive-tree, kd-tree, triangulation and spherical evaluation structures
  from `ev_atree.c`, `ev_kdtre.c`, `ev_trian.c`, `ev_sphere.c`, and the full
  evaluator orchestration in `ev_main.c`;
- process cross-validation from `procv.c`;
- simultaneous confidence-band/tube-formula machinery from `scb*.c` and
  `simul.c`;
- direct reproduction of the upstream internal QR/SVD/eigensystem APIs where
  the v0.1.0 public algorithms do not require them.

Those omissions are explicit so that v0.1.0 is not mistaken for a complete
replacement for every obscure `locfit` R feature.  No plotting code is planned
for translation.

## Numerical conventions

- Free-format source.
- `implicit none` throughout.
- Real kind is `dp = kind(1.0d0)`.
- Allocatable arrays replace upstream workspace pointer arithmetic.
- Explicit interfaces are supplied through modules.
- Dense linear algebra in the current core is self-contained and does not
  require BLAS/LAPACK.
