# Validation

The release tree was compiled with gfortran 14.2 using:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

and linked against LAPACK/BLAS.  All six permanent test programs and the example
pass.  FPM itself was not available in the validation container, so the exact
FPM source/dependency layout was compiled directly with the same compiler and
link libraries.

## Permanent tests

- `test_normal`: scalar/multivariate Gaussian and Student calculations,
  derivative tensors and simulation helpers.
- `test_binning`: generic linear binning, interpolation and convolution.
- `test_kde`: KDE/KDDE, CDF/quantile/random and grid helpers.
- `test_bandwidth_mixture`: normal-reference/LSCV/PI/SCV bandwidths and
  Gaussian-mixture MISE/AMISE/mode calculations.
- `test_highlevel`: KDA, 1D KCDE, deconvolution, histograms, curvature/modal
  significance and global two-sample testing.
- `test_boundary_support`: beta boundary KDE, empirical copula density,
  convex-hull/support logic, balloon KDE, and a multivariate KCDE reference.

## Independent randomized differential checks

### KDE / derivatives

120 randomized cases in dimensions 1--3 were compared with an independent
NumPy implementation using random SPD bandwidth matrices.  Maximum absolute
errors were:

- density: `1.6653345369377348e-16`
- gradient: `1.1102230246251565e-16`
- Hessian: `3.3306690738754696e-16`

### Multivariate KCDE

40 randomized two-dimensional data sets were compared with a direct SciPy
mixture of `multivariate_normal.cdf` evaluations.  With a tight probability
control, the maximum absolute discrepancy was `2.8236476843357394e-06` and the
mean absolute discrepancy was `1.8533342962574543e-07`.  The largest internal
reported Gaussian-CDF integration error was `4.02485840998215e-06`.

The permanent suite additionally contains a fixed 2D SciPy-reference KCDE case.

### Deconvolution weights

30 randomized small 2D problems were compared with an independent SciPy SLSQP
solution of the same simplex-constrained quadratic objective.  The Fortran
objective was never worse within recorded floating-point precision, the maximum
simplex-sum error was `2.22e-16`, and the maximum component-wise weight
difference was `8.41148e-06`.

## Release hygiene

The distributed tree contains no `.o`, `.mod`, executable, or temporary build
artifacts.  All translated Fortran is free-form `.f90` and uses explicit
interfaces for external numerical-library calls.
