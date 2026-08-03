! SPDX-License-Identifier: GPL-3.0-or-later
program test_pca
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rrcov, only : dp, pca_result, covariance_result, pca_classic, pca_cov, &
    pca_locantore, pca_grid, pca_hubert, cov_mcd, rrcov_success, rrcov_no_convergence
  implicit none
  real(dp) :: x(80, 4)
  type(pca_result) :: result
  type(covariance_result) :: robust_cov
  integer :: i

  do i = 1, 80
    x(i, 1) = sin(0.09_dp * real(i, dp))
    x(i, 2) = 2.0_dp * x(i, 1) + 0.10_dp * cos(0.31_dp * real(i, dp))
    x(i, 3) = -0.5_dp * x(i, 1) + 0.20_dp * sin(0.27_dp * real(i, dp))
    x(i, 4) = 0.15_dp * cos(0.41_dp * real(i, dp))
  end do
  x(1:4, :) = x(1:4, :) + 10.0_dp

  call pca_classic(x, result, k=2)
  call check_pca(result)
  call cov_mcd(x, robust_cov, nsamp=100, seed=10)
  call pca_cov(x, robust_cov, result, k=2)
  call check_pca(result)
  call pca_locantore(x, result, k=2)
  call check_pca(result)
  call pca_grid(x, result, k=2, max_directions=150, seed=20)
  call check_pca(result)
  call pca_hubert(x, result, k=2, kmax=3, nsamp=100, seed=30)
  call check_pca(result)

  print '(a)', "test_pca: PASS"
contains
  subroutine check_pca(value)
    type(pca_result), intent(in) :: value
    call assert_true(value%status == rrcov_success .or. value%status == rrcov_no_convergence, "PCA status")
    call assert_true(all(shape(value%loadings) == [4, 2]), "loading dimensions")
    call assert_true(all(shape(value%scores) == [80, 2]), "score dimensions")
    call assert_true(all(ieee_is_finite(value%loadings)), "finite loadings")
    call assert_true(all(ieee_is_finite(value%scores)), "finite scores")
    call assert_true(all(value%eigenvalues >= 0.0_dp), "nonnegative eigenvalues")
    call assert_true(size(value%score_distances) == 80, "score distances")
    call assert_true(size(value%orthogonal_distances) == 80, "orthogonal distances")
  end subroutine check_pca

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAIL: " // trim(message)
      error stop 1
    end if
  end subroutine assert_true
end program test_pca
