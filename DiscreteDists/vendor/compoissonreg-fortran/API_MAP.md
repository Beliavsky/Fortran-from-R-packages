# API map: COMPoissonReg 0.8.2 -> COMPoissonReg-fortran 0.1.0

## Exported R API

| R API | Fortran mapping | Notes |
|---|---|---|
| `dcmp`, `pcmp`, `qcmp`, `rcmp` | same base names | Scalar probability functions; RNG writes an integer array. |
| `ecmp`, `vcmp`, `ncmp`, `tcmp` | same names | Mean, variance, normalizer, truncation endpoint. |
| `dzicmp`, `pzicmp`, `qzicmp`, `rzicmp` | same base names | ZICMP probability/RNG functions. |
| `ezicmp`, `vzicmp` | same names | ZICMP moments. |
| `dzip`, `pzip`, `qzip`, `rzip` | same base names | ZIP probability/RNG functions. |
| `ezip`, `vzip` | same names | ZIP moments. |
| `glm.cmp.raw` | `fit_cmp_raw` | Matrix-based CMP MLE. |
| `glm.zicmp.raw` | `fit_zicmp_raw` | Matrix-based ZICMP MLE. |
| `glm.cmp` | not ported | R formula/model-frame parsing only; use the raw matrix interfaces. |
| `get.control` | `cmp_control_t` | Default type constructor gives upstream numerical defaults. |
| `get.init`, `get.init.zero` | `cmp_init_t`, `default_init` | Allocatable coefficient vectors. |
| `get.fixed` | `cmp_fixed_t`, `default_fixed` | Fortran uses logical fixed masks rather than R integer index vectors. |
| `get.offset`, `get.offset.zero` | `cmp_offset_t`, `default_offset` | Three offset vectors. |
| `get.modelmatrix` | direct matrix arguments | `X`, `S`, `W` are passed directly to fitting/prediction routines. |
| `equitest` | `equitest_cmp`, `equitest_zicmp` | Likelihood-ratio test. |
| `leverage` | `leverage_cmp` | The supplied R source contains no `leverage.zicmpfit` implementation. |
| `parametric.bootstrap` | `bootstrap_cmp`, `bootstrap_zicmp` | Returns coefficient replicates in caller storage. |
| `sdev` | `sdev_cmp`, `sdev_zicmp` | Standard errors from the inverse observed Hessian. |
| `nu` | `fitted_cmp` / `fitted_zicmp` | Deprecated upstream; link parameters are returned directly. |

## S3 computational methods

The R methods `AIC`, `BIC`, covariance, prediction, residuals, CMP leverage,
CMP deviance, and equidispersion testing map to `aic_*`, `bic_*`, the
`covariance` fit component, `predict_*`, `residuals_*`, `leverage_cmp`,
`deviance_cmp`, and `equitest_*`.

S3 printing and summary formatting are presentation code and are omitted.

## Native/internal computational routines

The following C++/R internals are also represented:

- truncated, asymptotic, and hybrid CMP normalizers;
- truncation endpoint;
- `z_prodj`, `z_prodj2`, `z_prodjlogj`, `z_prodlogj`, `z_prodlogj2`;
- CMP/ZICMP log-likelihoods;
- randomized quantile residual computation;
- `fim_cmp`, `fim_cmp_mc`, `fim_zicmp`, `fim_zicmp_mc`, `fim_zicmp_reg`.
