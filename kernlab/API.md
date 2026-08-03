# API overview

All public procedures are available through `use kernlab`.

## Kernel specifications

- `rbfdot(sigma)`
- `laplacedot(sigma)`
- `besseldot(sigma, order, degree)`
- `polydot(degree, scale, offset)`
- `tanhdot(scale, offset)`
- `vanilladot()`
- `anovadot(sigma, degree)`
- `splinedot()`
- `stringdot(length, lambda, normalized)`
- `fourierdot(sigma)`

These constructors return `type(kernel_spec)`.

## Kernel operations

- `kernel_value(kernel, x, y)`
- `kernel_matrix(kernel, x, K, status [, y])`
- `kernel_mult(kernel, x, z, result, status [, y, blocksize])`
- `kernel_pol(kernel, x, z, value, status [, y, k0])`
- `kernel_fast(kernel, x, y, values, status [, dota])`
- `string_kernel_value` and `string_kernel_matrix`

## Unsupervised and feature methods

- `kpca`, `kpca_predict`
- `kcca`
- `kha`, `kha_predict`
- `kfa`, `kfa_predict`
- `kkmeans`, `kkmeans_from_kernel`
- `specc`, `specc_from_kernel`
- `ranking`, `ranking_from_kernel`
- `csi`

## Supervised models

- Generic `ksvm` with integer targets for classification and real targets for
  regularized kernel regression.
- Generic `lssvm` with integer or real targets.
- Generic `gausspr` with integer or real targets.
- `rvm` for sparse Bayesian regression.
- `kqr` for kernel quantile regression.
- `onlearn` and `inlearn` for online kernel updates.
- `predict_kernel_model` for common prediction.
- `gausspr_predict_variance` for Gaussian-process posterior variances.

## Tests and utilities

- `kmmd`, `kmmd_from_kernels`
- `ipop`
- `inchol`
- `sigest`
- `couple`, `couple_vote`, `couple_pkpd`, `couple_minpair`

## Result types

- `kernel_spec`
- `kernel_model`
- `kpca_result`
- `kcca_result`
- `cluster_result`
- `mmd_result`
- `inchol_result`
- `ipop_result`
- `ranking_result`
- `csi_result`

Every non-scalar operation returns or sets an integer status. `KL_SUCCESS` is
zero; other status constants are declared in `kernlab_kinds` and re-exported by
`kernlab`.
