# Translation coverage

## Translated

- C++ `coneACpp` active-face polar-cone projection.
- C++ `coneBCpp` constraint-cone projection with an optional linear subspace.
- C++ `qprogCpp` Cholesky/cone-projection QP transformation.
- Weighted `coneA` and `coneB` transformations.
- Starting-face support.
- `qrdecomp` numerical column basis/rank.
- `makedelta` for all eight shape restrictions, including tied predictor values.
- `constreg` array-level constrained regression.
- The 2023 mixture covariance calculation in `constreg`.
- Coefficient/fit confidence intervals and t tests.
- Optional Monte-Carlo mixture-of-beta E01 tests.
- Array-level `shapereg.fit` for one constrained predictor plus optional parametric covariates.
- `check_irred` functionality as `check_irreducible`.
- Internal regularized-beta and Student-t CDF/quantile routines needed for inference.

## R-only code intentionally omitted

- Formula parsing and model-frame construction in `shapereg`.
- S3 `fitted`, `coef`, `summary`, `print`, and `vcov` dispatch/plumbing.
- Shape-marker functions whose only purpose is attaching R attributes (`incr`, `decr`, etc.); their eight numeric shape codes are public Fortran constants instead.
- Plotting and example graphics.
- R data objects and `.rda` serialization.

## Differences

1. The R/C++ implementation accepts Matrix-package sparse matrices but immediately converts them to dense Armadillo matrices. The Fortran API accepts dense arrays directly.
2. R validates weights as nonnegative, but its subsequent inverse-square-root transformation is singular for zero weights. The Fortran routines reject nonpositive weights explicitly with `coneproj_invalid_input`.
3. Armadillo solves are replaced by self-contained Cholesky/Gaussian/normal-equation solvers. A small ridge fallback is used only when an active-set Gram matrix is numerically singular.
4. The array-level `shapereg_fit` supports the computational `shapereg.fit` case of one shape-restricted predictor plus optional parametric columns. R formula/object orchestration for multiple decorated terms is not reproduced.
5. Monte-Carlo streams use the Fortran runtime RNG, so seeded simulated p-values/covariances are not bit-identical to R's RNG stream.
