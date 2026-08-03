! SPDX-License-Identifier: GPL-3.0-or-later
program test_core_estimators
  use intrinsicfrp, only: dp, vector_result, hj_result, screening_result
  use intrinsicfrp, only: tfrp, frp, sdf_coefficients, hj_misspecification_distance
  use intrinsicfrp, only: gkr_factor_screening
  use intrinsicfrp, only: covariance_matrix, cross_covariance, hac_covariance
  implicit none
  integer, parameter :: n = 240, p = 6, k = 3
  real(dp) :: factors(n, k), returns(n, p), beta(p, k), alpha(p)
  real(dp), allocatable :: cov_fr(:, :), var_r(:, :), mean_r(:), manual(:)
  real(dp), allocatable :: temp(:, :), rhs(:, :), hac(:, :)
  type(vector_result) :: out
  type(hj_result) :: hj
  type(screening_result) :: screen
  integer :: i, j, st

  do i = 1, n
    factors(i, 1) = sin(0.071_dp * real(i, dp))
    factors(i, 2) = cos(0.043_dp * real(i, dp))
    factors(i, 3) = sin(0.019_dp * real(i, dp)) + 0.2_dp * cos(0.11_dp * real(i, dp))
  end do
  do j = 1, p
    alpha(j) = 0.002_dp * real(j, dp)
    beta(j, 1) = 0.25_dp + 0.03_dp * real(j, dp)
    beta(j, 2) = -0.15_dp + 0.02_dp * real(j, dp)
    beta(j, 3) = 0.08_dp * real(mod(j, 3) - 1, dp)
  end do
  do i = 1, n
    returns(i, :) = alpha + matmul(beta, factors(i, :)) + &
      0.025_dp * [(sin(0.17_dp * real(i + 3 * j, dp)), j = 1, p)]
  end do

  call tfrp(returns, factors, out, include_standard_errors=.true.)
  call check(out%status == 0, 'tfrp status')
  call check(size(out%estimate) == k, 'tfrp length')
  call check(size(out%standard_errors) == k, 'tfrp se length')
  call check(all(out%standard_errors >= 0.0_dp), 'tfrp se nonnegative')

  cov_fr = cross_covariance(factors, returns)
  var_r = covariance_matrix(returns)
  allocate(mean_r(p), rhs(p, 1))
  mean_r = sum(returns, dim=1) / real(n, dp)
  rhs(:, 1) = mean_r
  call solve_local(var_r, rhs, temp, st)
  manual = matmul(cov_fr, temp(:, 1))
  call check(maxval(abs(out%estimate - manual)) < 1.0e-9_dp, 'tfrp formula')

  call frp(returns, factors, out, misspecification_robust=.false., &
    include_standard_errors=.true.)
  call check(out%status == 0 .and. size(out%estimate) == k, 'fm frp')
  call check(all(abs(out%estimate) < 10.0_dp), 'fm frp finite scale')

  call frp(returns, factors, out, misspecification_robust=.true., &
    include_standard_errors=.true.)
  call check(out%status == 0 .and. all(out%standard_errors >= 0.0_dp), 'krs frp')

  call sdf_coefficients(returns, factors, out, misspecification_robust=.false., &
    include_standard_errors=.true.)
  call check(out%status == 0 .and. size(out%estimate) == k, 'fm sdf')

  call sdf_coefficients(returns, factors, out, misspecification_robust=.true., &
    include_standard_errors=.true.)
  call check(out%status == 0 .and. size(out%estimate) == k, 'gkr sdf')


  call gkr_factor_screening(returns, factors, screen, target_level=0.05_dp)
  call check(screen%status == 0, 'gkr screening status')
  call check(size(screen%selected_indices) <= k, 'gkr screening indices')

  call frp(returns, factors, out, misspecification_robust=.true., &
    screening_level=0.05_dp)
  call check(out%status == 0, 'frp screening')

  call hj_misspecification_distance(returns, factors, hj)
  call check(hj%status == 0, 'hj status')
  call check(hj%squared_distance >= -1.0e-8_dp, 'hj nonnegative')
  call check(hj%lower_bound <= hj%upper_bound, 'hj interval order')

  call hac_covariance(returns - spread(mean_r, 1, n), hac, st, prewhite=.true.)
  call check(st == 0 .and. size(hac, 1) == p, 'hac dimensions')
  call check(maxval(abs(hac - transpose(hac))) < 1.0e-10_dp, 'hac symmetric')

  print '(a)', 'test_core_estimators: PASS'

contains

  subroutine solve_local(a, b, x, status)
    use intrinsicfrp_linalg, only: solve_linear
    real(dp), intent(in) :: a(:, :), b(:, :)
    real(dp), allocatable, intent(out) :: x(:, :)
    integer, intent(out) :: status
    call solve_linear(a, b, x, status)
  end subroutine solve_local

  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a,1x,a)') 'FAIL:', trim(label)
      error stop 1
    end if
  end subroutine check
end program test_core_estimators
