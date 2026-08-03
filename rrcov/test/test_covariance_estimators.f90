! SPDX-License-Identifier: GPL-3.0-or-later
program test_covariance_estimators
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rrcov, only : dp, covariance_result, cov_classic, cov_mcd, cov_mve, cov_ogk, &
    cov_mest, cov_sest, cov_mmest, cov_sde, cov_mrcd, rrcov_success, rrcov_no_convergence
  implicit none
  real(dp) :: x(60, 3), classic_norm, robust_norm
  type(covariance_result) :: classic, mcd, estimate
  integer :: i

  do i = 1, 60
    x(i, 1) = sin(0.17_dp * real(i, dp)) + 0.10_dp * cos(0.43_dp * real(i, dp))
    x(i, 2) = 0.7_dp * x(i, 1) + 0.4_dp * cos(0.11_dp * real(i, dp))
    x(i, 3) = -0.3_dp * x(i, 1) + 0.5_dp * sin(0.29_dp * real(i, dp))
  end do
  x(1:6, 1) = x(1:6, 1) + 15.0_dp
  x(1:6, 2) = x(1:6, 2) - 12.0_dp
  x(1:6, 3) = x(1:6, 3) + 8.0_dp

  call cov_classic(x, classic)
  call check_estimate(classic)
  call cov_mcd(x, mcd, nsamp=120, seed=1234)
  call check_estimate(mcd)
  classic_norm = sqrt(sum(classic%center ** 2))
  robust_norm = sqrt(sum(mcd%center ** 2))
  call assert_true(robust_norm < classic_norm, "MCD center should resist gross outliers")

  call cov_mve(x, estimate, nsamp=120, seed=4321)
  call check_estimate(estimate)
  call cov_ogk(x, estimate, niter=2)
  call check_estimate(estimate)
  call cov_mest(x, estimate)
  call check_estimate(estimate)
  call cov_sest(x, estimate, nsamp=80, seed=99)
  call check_estimate(estimate)
  call cov_mmest(x, estimate, nsamp=80, seed=99)
  call check_estimate(estimate)
  call cov_sde(x, estimate, ndirections=120, seed=77)
  call check_estimate(estimate)
  call cov_mrcd(x, estimate, nsamp=80, seed=55)
  call check_estimate(estimate)

  print '(a)', "test_covariance_estimators: PASS"
contains
  subroutine check_estimate(value)
    type(covariance_result), intent(in) :: value
    call assert_true(value%status == rrcov_success .or. value%status == rrcov_no_convergence, &
      "estimator status")
    call assert_true(size(value%center) == 3, "center dimension")
    call assert_true(all(shape(value%covariance) == [3, 3]), "covariance dimension")
    call assert_true(all(ieee_is_finite(value%center)), "finite center")
    call assert_true(all(ieee_is_finite(value%covariance)), "finite covariance")
    call assert_true(maxval(abs(value%covariance - transpose(value%covariance))) < 1.0e-8_dp, &
      "symmetric covariance")
    call assert_true(all([(value%covariance(i, i) > 0.0_dp, i=1, 3)]), "positive diagonal")
  end subroutine check_estimate

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAIL: " // trim(message)
      error stop 1
    end if
  end subroutine assert_true
end program test_covariance_estimators
