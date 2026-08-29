# API mapping

| R / upstream operation | Fortran API / implementation |
|---|---|
| spam/as.spam/triplet/as.matrix | `csr_matrix`, `csr_from_dense`, `csr_from_triplet`, `csr_to_dense` |
| diag.spam / diag<- | `csr_diag`, `csr_diagonal`, `csr_set_diagonal` |
| t.spam | `csr_transpose` |
| +, -, element multiplication | `csr_add`, `csr_subtract`, `csr_hadamard` |
| matrix multiplication | `csr_matmul`, `csr_matvec`, `csr_matmat` |
| kronecker.spam | `csr_kronecker` |
| rbind.spam/cbind.spam | `csr_rbind`, `csr_cbind` |
| subset.spam | `csr_subset` |
| upper.tri/lower.tri structural extraction | `csr_upper`, `csr_lower` |
| rowSums/colSums/rowMeans/colMeans | `csr_row_sums`, `csr_col_sums`, `csr_row_means`, `csr_col_means` |
| norm.spam | `csr_norm` |
| bandwidth | `csr_bandwidth` |
| var.spam | `csr_var` |
| diff.spam | `csr_diff` |
| toeplitz.spam/circulant.spam | `csr_toeplitz`, `csr_circulant` |
| bdiag.spam | `csr_bdiag` |
| crossprod/tcrossprod | `csr_crossprod`, `csr_tcrossprod` |
| apply.spam | `csr_map_entries`, `csr_apply_margin` typed callbacks |
| permutation.spam | `permute_csr`, `inverse_permutation` |
| gmult | `gmult` |
| spam_random | `spam_random_matrix` |
| nearest.dist | `nearest_dist` |
| spam_rdist | `rdist` |
| spam_rdist.earth | `rdist_earth` |
| cov.exp | `cov_exp` |
| cov.sph/cor.sph | `cov_sph`, `cor_sph` |
| cov.nug | `cov_nug` |
| cov.wu1/2/3 | `cov_wu1`, `cov_wu2`, `cov_wu3` |
| cov.wend1/2 | `cov_wend1`, `cov_wend2` |
| cov.mat / cov.finnmat | `cov_mat`, `cov_finnmat` |
| cov.mat12/32/52 | `cov_mat12`, `cov_mat32`, `cov_mat52` |
| precmat.RW1/RW2/RWn | `precmat_rw1`, `precmat_rw2`, `precmat_rwn` |
| precmat.season | `precmat_season` |
| precmat.IGMRFreglat/IGMRFirreglat | `precmat_igmrf_reglat`, `precmat_igmrf_irreglat` |
| precmat.GMRFreglat | `precmat_gmrf_reglat` |
| chol.spam | `spam_chol_factor` -> `spam_chol` |
| update.spam.chol.NgPeyton | `spam_chol_update` |
| solve.spam | `spam_solve` |
| forwardsolve/backsolve | `spam_forwardsolve`, `spam_backsolve` |
| determinant.spam | `spam_logdet`, `spam_determinant` |
| determinant(chol factor) | `chol_factor_logdet` |
| chol2inv.spam | `spam_chol2inv` |
| eigen.spam/eigen_approx | `spam_eigen_symmetric`, `spam_eigen_general` (bundled ARPACK) |
| rmvnorm(.spam) | `rmvnorm_cov` |
| rmvnorm.prec | `rmvnorm_prec` |
| rmvnorm.canonical | `rmvnorm_canonical` |
| constrained Normal variants | `rmvnorm_cov_const`, `rmvnorm_prec_const`, `rmvnorm_canonical_const` |
| rmvnorm.conditional | `rmvnorm_conditional` |
| rmvt(.spam) | `rmvt_cov` |
| rgrf | `rgrf`, `rgrf_grid` |
| neg2loglikelihood.spam | `neg2loglikelihood_spam` |
| neg2loglikelihood.nomean | `neg2loglikelihood_nomean_spam` |
| mle.spam/mle | `mle_spam` for built-ins; `mle_spam_custom` for arbitrary covariance callbacks |
| mle.nomean(.spam) | `mle_nomean_spam`; `mle_nomean_spam_custom` |
| SparseKit native routines | retained under `src/native/fromsparsekit.f90` and related translated sources |
| Ng-Peyton/ordering kernels | retained under `src/native/cholmodified.f90`, `bckslvmodified.f90`, `permutation.f90` |
| ARPACK | retained under `src/native/ds_ARPACK.f90`, `dn_ARPACK.f90`, `dgetv0.f90` |

R S4/S3 dispatch, plotting/display, Matrix/dotCall64/Rcpp interop, package option helpers, and foreign sparse file readers are intentionally omitted as non-computational/interface code.
