# Porting notes

## Upstream

- Package: DirichletReg
- Upstream version: 0.7-2
- Date in DESCRIPTION: 2025-05-30
- Author: Marco Johannes Maier
- Upstream license: GPL (>= 2)

The original `DESCRIPTION`, `NAMESPACE`, `NEWS`, R sources, and C sources are
retained in `upstream/` for attribution and auditability.

## Computational mapping

| Upstream code | Fortran implementation |
| --- | --- |
| `R/Dirichlet.R`, `src/ddirichlet.c`, `src/rdirichlet.c` | `dirichletreg_distribution.f90`, `dirichletreg_rng.f90` |
| `R/DR_data.R` | `dirichletreg_data.f90` |
| `R/DR_LL_common.R`, `src/wght_LL_grad_common.c` | `dirichletreg_common.f90` |
| `R/DR_LL_alt.R`, `src/wght_LL_grad_alternative.c` | `dirichletreg_alternative.f90` |
| `R/get_starting_values.R` | `dirichletreg_start.f90` |
| `R/DirichReg_fit.R` | `dirichletreg_fit.f90`, `dirichletreg_optimize.f90` |
| `fitted` / `predict` methods | `common_predict`, `alternative_predict`, fitted arrays in `dirichletreg_model` |
| `residuals.DirichletRegModel.R` | `dirichletreg_inference.f90` |
| `confint`, AIC, BIC, ANOVA numerics | model fields plus `dirichletreg_inference.f90` |
| `toTernaryQuaternary.R` | `dirichletreg_geometry.f90` |

## Omitted R-only or plotting functionality

- All 2-D/3-D/4-D plotting and color-conversion helpers.
- R formula parsing, `Formula` dependency, model frames, S3 print/summary
  presentation, `terms`, `model.matrix`, and `update` machinery.
- Formula-driven `drop1`. The lower-level likelihood-ratio calculation is
  implemented.
- Bundled R datasets are not converted to a Fortran data serialization format.

## Intentional corrections / numerical changes

### Starting-value beta score

The upstream `get_starting_values.R` beta-regression derivative uses
`psigamma(exp(eta1 + eta2))` and omits the chain-rule alpha/beta multiplier.
For a beta model with `a = exp(eta1)` and `b = exp(eta2)`, the corresponding
score uses `digamma(a+b)` and includes the multiplier `a` or `b`. The Fortran
port uses that mathematically consistent score. This routine is only used to
obtain starting values; the fitted Dirichlet target likelihood is unchanged.

### Unused `src/wght_LL.c`

The upstream tree contains an unregistered and apparently unused `wght_LL.c`.
Its scalar accumulator is not reset per observation before applying a weight,
which would make it unsuitable as a direct weighted likelihood implementation.
The active upstream code uses the registered gradient routines instead. The
Fortran port therefore implements the correct weighted likelihood directly in
both active regression evaluators rather than reproducing this unused routine.

### Stable softmax

The alternative model computes the same softmax means as the R code but
subtracts the row maximum before exponentiation. This is algebraically
identical and avoids overflow for large linear predictors.

### Optimizer

The R package uses `maxLik::maxBFGS` followed by `maxLik::maxNR`. The Fortran
port implements the same BFGS -> Newton strategy internally, with backtracking
line searches. Therefore iteration counts and the last few floating-point bits
can differ even when the optimum is the same.

## Validation

The Fortran tests include:

- independent Dirichlet log-density reference values;
- RNG normalization and empirical-mean checks;
- composition normalization / boundary transformation checks;
- central finite-difference validation of the common analytic score/Hessian;
- central finite-difference validation of the alternative analytic
  score/Hessian;
- simulated MLE recovery for common and alternative intercept-only models.

These are designed to test the numerical kernels independently of R's formula
and object system.
