module irlba
  use irlba_kinds, only : dp, i4
  use irlba_sparse, only : csc_matrix, csc_from_dense
  use irlba_core, only : irlba_result, irlba_control, irlb_operator
  use irlba_operator, only : linear_operator, dense_operator, csc_operator
  use irlba_algorithms, only : irlba_svd => irlba, irlba_dense, irlba_sparse_matrix, irlba_complex, &
    partial_eigen, prcomp_irlba, ssvd, svdr, complex_svd_result, eigen_result, pca_result, ssvd_result, svdr_result
  implicit none
  public
end module irlba
