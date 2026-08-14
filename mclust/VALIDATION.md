# Validation

Validation was performed with GNU Fortran 14.2.0 and external BLAS/LAPACK.
Modern free-form modules and tests were compiled with:

```text
-std=f2018 -O2 -fcheck=all -Wall -Wextra -Werror
-Wno-maybe-uninitialized -Wimplicit-interface
```

The upstream numerical kernels, converted from fixed source form to free source form without algorithmic modernization, were compiled separately
with `-O2`; the modern callers use explicit interfaces and runtime checks.

## Permanent tests

Seven permanent test programs pass:

1. `test_models` -- EM over the full model family and BIC selection.
2. `test_covariance_constraints` -- checks the defining volume/shape/orientation
   invariants of all 14 multivariate covariance models.
3. `test_density_sim` -- density, one-dimensional CDF/quantiles, HDR, and
   simulation.
4. `test_classification_ssc` -- discriminant analysis and semi-supervised
   mixture fitting.
5. `test_utilities` -- imputation, MclustDR, weighted fitting, cluster
   combination, and numerical utilities.
6. `test_bootstrap` -- parametric bootstrap paths.
7. `test_mapping_crimcoords` -- unmap, label matching, majority voting, and
   canonical discriminant coordinates.

The example also passes and selects a two-component EEE model on its synthetic
data.

## Full-pipeline diabetes reference

The translated default model-selection pipeline was run on the 145-observation
mclust diabetes example after applying the data correction shown in the
upstream rendered vignette.

Upstream vignette result:

- model: VVV
- components: 3
- cluster sizes: 82, 35, 28
- mixing proportions: 0.5553630, 0.2479432, 0.1966939
- log likelihood: -2295.118
- BIC: -4734.561
- ICL: -4749.16

Fortran result:

- model: VVV
- components: 3
- cluster sizes: 82, 35, 28
- mixing proportions: 0.5553629738, 0.2479431588, 0.1966938674
- log likelihood: -2295.10858565
- BIC: -4734.54244983
- ICL: -4749.04017667

Thus model, component count, classifications, and displayed mixing proportions
agree; the objective values differ slightly because the modern high-level EM
loop has its own convergence path around the retained covariance M-steps.
`validation/diabetes_reference.f90` turns these values into regression checks.

## Independent Gaussian-density validation

`dmvnorm` was compared against SciPy's multivariate-normal log density on 250
random positive-definite covariance matrices in dimensions 1 through 6.  The
maximum absolute log-density difference was approximately **4.51e-13**.

## Build hygiene

- `fpm.toml` parses as TOML.
- Source/test/example Fortran lines are at most 132 columns, so default free-form
  line-length rules do not require a compiler-specific override.
- The final archive is generated from the source tree only; `.o`, `.mod`,
  executables, caches, and temporary validation files are excluded.

Run `validation/run_strict.sh` on a Unix-like system with gfortran, BLAS, and
LAPACK to reproduce the strict build and permanent tests.

## FPM compatibility fix in v0.1.2

The project manifest now sets `implicit-typing = true`, `implicit-external = true`, and `source-form = "free"`. This is required because fpm disables implicit typing and implicit external interfaces by default, while the translated historical MCLUST kernels intentionally retain those valid legacy language features. The full regression suite was rerun after this manifest change with no numerical changes.
