# Dependencies

The target layout is the root of `Fortran-from-R-packages`, with these sibling directories:

```text
Fortran-from-R-packages/
  kSamples/
  rfortran-core/
  SuppDists/
```

`fpm.toml` uses:

```toml
[dependencies]
rfortran-core = { path = "../rfortran-core" }
suppdists-fortran = { path = "../SuppDists" }
```

`rfortran-core` supplies `r_kinds::dp`, normal and chi-square distribution functions, average ranks, and sample variance. `SuppDists` supplies `suppdists::norm_order` for the normal-score QN variant.

The translated package has no direct BLAS/LAPACK/ARPACK requirement and therefore does not link system numerical libraries or copy any such implementation.
