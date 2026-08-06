# API map

| Upstream R routine | Fortran routine | Coverage |
|---|---|---|
| `StLM` | `stlm`, `fit_student_lm` | Static Student-t regression |
| `StAR` | `star`, `fit_student_ar` | Univariate autoregression |
| `StDLM` | `stdlm`, `fit_student_dlm` | Dynamic linear regression |
| `StVAR` | `stvar`, `fit_student_var` | Vector autoregression |
| `BlockTop` | `block_top_covariance` | Exact block-factor construction |
| `Par.dlrm`, `Par.star`, `Par.stvar` | `conditional_parameters`, internal decoder | Conditional parameters |
| `Jacob.*` | internal central-difference Jacobian | Coefficient inference |
| `ConJacob.*` | internal central-difference Jacobian | Variance-coefficient inference |
| `Student` | internal Student-t AD diagnostic | Standardization and AD test |
| `ADGofTest::ad.test` | self-contained Marsaglia AD approximation | No external dependency |
| `numDeriv::jacobian` | self-contained central differences | No external dependency |
| `MCMCpack::vech`, `xpnd` | `vech`, `expand_vech` | Lower-triangle conversion |
| `matlab::ones` | intrinsic array assignment | No dependency required |
| `stats::optim` | self-contained BFGS/Nelder-Mead | Numerical minimization |
| `tseries::get.hist.quote` | omitted | Imported but unused upstream |
