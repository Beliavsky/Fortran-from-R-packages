! SPDX-License-Identifier: GPL-3.0-or-later
program demo_rrcov
  use rrcov, only : dp, covariance_result, pca_result, lda_model, &
    cov_classic, cov_mcd, pca_cov, lda_cov_fit, lda_predict, outlier_flags
  implicit none
  real(dp) :: x(100, 3)
  integer :: grouping(100), i
  integer, allocatable :: predicted(:)
  logical, allocatable :: flags(:)
  type(covariance_result) :: classic, robust
  type(pca_result) :: pca
  type(lda_model) :: lda

  do i = 1, 50
    grouping(i) = 1
    x(i, 1) = -1.5_dp + 0.5_dp * sin(0.17_dp * real(i, dp))
    x(i, 2) = -0.8_dp + 0.4_dp * cos(0.23_dp * real(i, dp))
    x(i, 3) = 0.3_dp * x(i, 1) + 0.2_dp * sin(0.31_dp * real(i, dp))
  end do
  do i = 51, 100
    grouping(i) = 2
    x(i, 1) = 1.5_dp + 0.5_dp * sin(0.19_dp * real(i, dp))
    x(i, 2) = 0.8_dp + 0.4_dp * cos(0.29_dp * real(i, dp))
    x(i, 3) = 0.3_dp * x(i, 1) + 0.2_dp * sin(0.37_dp * real(i, dp))
  end do
  x(1:5, :) = x(1:5, :) + 10.0_dp

  call cov_classic(x, classic)
  call cov_mcd(x, robust, nsamp=250, seed=2026)
  flags = outlier_flags(robust)
  call pca_cov(x, robust, pca, k=2)
  call lda_cov_fit(x, grouping, lda, "mcd", nsamp=150, seed=2026)
  call lda_predict(lda, x, predicted)

  print '(a,3f11.5)', "Classical center: ", classic%center
  print '(a,3f11.5)', "Robust MCD center:", robust%center
  print '(a,i0)', "Flagged robust outliers: ", count(flags)
  print '(a,2f11.5)', "Robust PCA eigenvalues: ", pca%eigenvalues
  print '(a,f8.3,a)', "Training classification accuracy: ", &
    100.0_dp * real(count(predicted == grouping), dp) / real(size(grouping), dp), "%"
end program demo_rrcov
