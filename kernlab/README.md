# kernlab-fortran

A self-contained modern Fortran/FPM translation of the computational surface of
R package **kernlab 0.9-33**.

## Implemented areas

- Numeric kernels: Gaussian RBF, Laplace, Bessel, polynomial, hyperbolic
  tangent, linear, ANOVA, spline, and Fourier.
- Normalized fixed-length spectrum string kernel.
- Kernel matrices, matrix products, quadratic forms, and fast kernel columns.
- Kernel PCA, kernel CCA, kernel Hebbian/batch subspace extraction, and kernel
  feature-basis selection.
- Kernel k-means and spectral clustering.
- SVM classification, regularized kernel regression, least-squares SVM,
  Gaussian-process regression, relevance-vector regression, and kernel
  quantile regression.
- Manifold ranking, online kernel learning, kernel MMD tests, incomplete
  Cholesky, CSI reduction, RBF bandwidth estimation, probability coupling, and
  a box/two-sided constrained quadratic-programming interface.

R formula methods, S4 classes, plotting, display methods, bundled data, and
R-specific prediction dispatch are not compiled. They are replaced by explicit
array interfaces and derived result types.

## Build

```text
fpm build
fpm test
fpm run demo_kernlab
```

The library uses `real(dp)`, where `dp = kind(1.0d0)`, and has no external
Fortran package dependencies.

## Minimal example

```fortran
use kernlab
real(dp) :: x(8,2)
integer :: y(8), status
type(kernel_spec) :: kernel
type(kernel_model) :: model
real(dp), allocatable :: scores(:,:)
integer, allocatable :: predicted(:)

kernel = rbfdot(0.5_dp)
call ksvm(x, y, kernel, model, cost=10.0_dp)
call predict_kernel_model(model, x, scores, status, predicted)
```

See `API.md`, `PORTING.md`, and the programs under `example/`.
