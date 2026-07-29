! SPDX-License-Identifier: GPL-3.0-or-later
program test_risk
  use frapo
  implicit none

  real(dp) :: covariance(2, 2), weights(2), expected_sigma
  real(dp) :: x(6, 3)
  real(dp), allocatable :: contribution(:), dependence(:, :)
  integer :: status

  covariance = reshape([4.0_dp, 1.0_dp, 1.0_dp, 9.0_dp], [2, 2])
  weights = [0.6_dp, 0.4_dp]
  expected_sigma = sqrt(dot_product(weights, matmul(covariance, weights)))

  call assert_close(dr(weights, covariance), 2.4_dp / expected_sigma, 1.0e-13_dp, 'DR')
  call assert_close(cr(weights, covariance), 0.5_dp, 1.0e-14_dp, 'CR')
  call assert_close(rhow(weights, covariance), 1.0_dp / 6.0_dp, 1.0e-14_dp, 'rho-w')

  contribution = mrc(weights, covariance, percentage=.false., status=status)
  call assert_equal_int(status, frapo_ok, 'MRC status')
  call assert_close(sum(contribution), expected_sigma, 1.0e-13_dp, 'MRC decomposition')
  contribution = mrc(weights, covariance, percentage=.true.)
  call assert_close(sum(contribution), 100.0_dp, 1.0e-12_dp, 'MRC percentages')

  x(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
  x(:, 2) = x(:, 1)
  x(:, 3) = [6.0_dp, 5.0_dp, 4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]
  dependence = tdc(x, method=tdc_empirical, lower_tail=.true., k=2, status=status)
  call assert_equal_int(status, frapo_ok, 'TDC status')
  call assert_close(dependence(1, 2), 1.0_dp, 0.0_dp, 'concordant empirical TDC')
  call assert_close(dependence(1, 3), 0.0_dp, 0.0_dp, 'opposite empirical TDC')
  dependence = tdc(x, method=tdc_evt, lower_tail=.false., k=2)
  call assert_close(dependence(1, 2), 1.0_dp, 0.0_dp, 'concordant EVT TDC')
  call assert_close(dependence(1, 3), 0.0_dp, 0.0_dp, 'opposite EVT TDC')

  print '(a)', 'test_risk: PASS'

contains
  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tolerance) then
      write(*, '(a,2es24.15)') trim(label)//': ', actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_equal_int(actual, expected, label)
    integer, intent(in) :: actual, expected
    character(len=*), intent(in) :: label
    if (actual /= expected) then
      write(*, '(a,2i8)') trim(label)//': ', actual, expected
      error stop 1
    end if
  end subroutine assert_equal_int
end program test_risk
