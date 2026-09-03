# API coverage

Source snapshot: R package **fda 6.3.0**, dated 2025-05-21.

This document distinguishes translated numerical algorithms from R object,
dispatch, plotting, formula, and example/data-set infrastructure.

## Basis and penalty layer

| Upstream R API | Fortran API | Coverage |
|---|---|---|
| `create.constant.basis` | `make_constant_basis` | Computational equivalent |
| `create.bspline.basis`, `bsplineS` | `make_bspline_basis`, `eval_basis` | B-spline construction/evaluation and integer derivatives |
| `create.fourier.basis`, `fourier` | `make_fourier_basis`, `eval_basis` | Normalized Fourier values and derivatives |
| `create.monomial.basis`, `monomial` | `make_monomial_basis`, `eval_basis` | Integer powers; mathematically correct derivative scaling |
| `create.exponential.basis`, `expon` | `make_exponential_basis`, `eval_basis` | Values and integer derivatives |
| `create.power.basis`, `powerbasis` | `make_power_basis`, `eval_basis` | Real powers and integer derivatives where defined |
| `create.polygonal.basis`, `polyg` | `make_polygonal_basis`, `eval_basis` | Piecewise-linear values and first derivatives |
| `eval.basis`, `getbasismatrix` | `eval_basis` | Basis evaluation for the seven types above |
| `eval.penalty`, `*pen` routines | `basis_penalty`, `basis_gram` | Integer-derivative L2 penalties by composite Simpson quadrature |

Not yet covered: arbitrary `Lfd` coefficient-function differential operators,
`dropind`, cached R basis-value lists, and every specialized analytic penalty
shortcut.  The general numerical penalty produces the same mathematical
quantity for the translated integer-derivative case.

## Functional data and integration

| Upstream R API | Fortran API | Coverage |
|---|---|---|
| `fd` | `make_fd`, `fd_type` | Univariate coefficient matrices |
| `eval.fd`, `predict.fd` numerical core | `eval_fd` | Values and integer derivatives |
| `mean.fd` | `mean_fd` | Coefficientwise replication mean |
| `center.fd` | `center_fd` | Replication centering |
| `inprod` | `inprod_basis`, `inprod_fd` | Cross-basis and cross-replication L2 products |
| `trapzmat` | `trapz_mat` | Weighted equal-grid trapezoidal matrix products |
| `quadset` | `quadset` | Composite Simpson points and weights |
| `polintmat` | `polint_matrix` | Neville polynomial interpolation/extrapolation for matrix-valued sequences |
| `symsolve` | `symsolve` | Symmetric positive-definite solves |
| `zerofind` | `zero_find` | Range-bracketing zero test |

Not yet covered: three-dimensional/multivariate `fd` arrays, `bifd`, R names and
class methods, overloaded arithmetic/object concatenation methods, and `lnsrch`
(the callback-based R helper is omitted rather than violating the project rule
that every Fortran data dummy carry `INTENT` or `VALUE`).

## Smoothing and model-selection layer

| Upstream R API | Fortran API | Coverage |
|---|---|---|
| `smooth.basis`, `smooth.basis1` | `smooth_basis` | Dense univariate matrix `y`, vector weights, roughness penalty, SSE/df/GCV/map |
| `lambda2df` | `lambda_to_df` | Effective df from the smoother trace |
| `df2lambda` | `df_to_lambda` | Monotone log-scale inversion by bisection |
| `lambda2gcv` | `lambda_to_gcv` | Per-replication GCV |
| `project.basis` | `project_basis` | Least-squares projection plus optional small stabilizing penalty |

Not yet covered: matrix-valued observation weights, semiparametric covariate
augmentation, three-dimensional response arrays, the upstream QR branch, and
full `fdPar`/`Lfd` object semantics.

## Functional PCA and CCA

| Upstream R API | Fortran API | Coverage |
|---|---|---|
| `pca.fd` | `pca_fd`, `pca_result_type` | Univariate regularized FPCA with optional centering |
| `cca.fd` | `cca_fd`, `cca_result_type` | Univariate regularized functional CCA |
| `geigen` | `geigen` | Generalized SVD with SPD row/column metrics |

The translated routines accept explicit bases, derivative orders, and lambdas
rather than reproducing R `fdPar` objects.

## ODE layer

| Upstream R API | Fortran API | Coverage |
|---|---|---|
| `derivs` | `linear_ode_rhs` | Order-`m` homogeneous operator to first-order system |
| `odesolv` | `odesolv`, `ode_solution_type` | Adaptive Cash–Karp integration of one or more initial conditions |
| `rkqs`, `rkck` | private Fortran routines | Adaptive embedded-step machinery |

Intentional correction: upstream `rkck` contains `575/512` for the second stage
in its sixth-stage state.  The translated implementation uses the standard
Cash–Karp coefficient `175/512`; deterministic testing against `u''+u=0`
confirms the corrected solver reaches the expected endpoint accurately.

`monomial` is likewise translated according to its documented mathematical
specification rather than two apparent defects in the supplied R source: the
upstream derivative recurrence has an off-by-one falling-factorial term and
omits the derivative scale induced by `argtrans`.  The Fortran version uses the
correct falling factorial and argument scaling; this is recorded in
`NOTICE.md`.

## Major meaningful parity targets not yet translated

- registration/warping and amplitude-phase algorithms (`register.fd`,
  landmark registration, monotone warps, phase-amplitude decomposition);
- PACE/sparse-longitudinal covariance and FPCA workflows;
- functional linear regression, regression inference, `linmod`, and related
  model-matrix/formula orchestration;
- principal differential analysis (`pda`, `eigen.pda`) and several specialized
  differential-operator estimation workflows;
- bivariate/multivariate functional-data (`bifd`) object algebra;
- specialized smoothing families such as monotone/intensity/positive smooths;
- R S3 methods, printing/summaries, plotting, formula evaluation, data loading,
  vignettes, and example-specific presentation code.

These omissions are explicit parity gaps, not stub implementations.
