! SPDX-License-Identifier: GPL-3.0-or-later
program test_statistics
  use rrcov, only : dp, median, mad_scale, qn_scale, sn_scale, tau_scale, &
    chi_square_cdf, chi_square_quantile, covariance_matrix, matrix_sqrt, &
    covariance_to_correlation, ilr_transform, rrcov_success
  implicit none
  real(dp) :: x(5), data(5, 2), probability, q
  real(dp), allocatable :: covariance(:, :), root(:, :), correlation(:, :), composition(:, :), transformed(:, :)
  integer :: status
  x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 100.0_dp]
  call assert_close(median(x), 3.0_dp, 1.0e-12_dp, "median")
  call assert_true(mad_scale(x) > 0.0_dp, "MAD")
  call assert_true(qn_scale(x) > 0.0_dp, "Qn")
  call assert_true(sn_scale(x) > 0.0_dp, "Sn")
  call assert_true(tau_scale(x) > 0.0_dp, "tau")
  q = chi_square_quantile(0.95_dp, 3.0_dp)
  probability = chi_square_cdf(q, 3.0_dp)
  call assert_close(probability, 0.95_dp, 1.0e-8_dp, "chi-square inverse")

  data(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
  data(:, 2) = [2.0_dp, 1.0_dp, 4.0_dp, 3.0_dp, 6.0_dp]
  covariance = covariance_matrix(data, status=status)
  root = matrix_sqrt(covariance, status=status)
  call assert_true(maxval(abs(matmul(root, root) - covariance)) < 1.0e-7_dp, "matrix square root")
  correlation = covariance_to_correlation(covariance, status)
  call assert_close(correlation(1, 1), 1.0_dp, 1.0e-12_dp, "correlation diagonal")

  allocate(composition(2, 3))
  composition(1, :) = [0.2_dp, 0.3_dp, 0.5_dp]
  composition(2, :) = [0.4_dp, 0.4_dp, 0.2_dp]
  call ilr_transform(composition, transformed, status)
  call assert_true(status == rrcov_success, "ILR status")
  call assert_true(all(shape(transformed) == [2, 2]), "ILR dimension")

  print '(a)', "test_statistics: PASS"
contains
  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    if (abs(actual - expected) > tolerance) then
      write(*, '(a,2es16.7)') "FAIL: " // trim(message) // " ", actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAIL: " // trim(message)
      error stop 1
    end if
  end subroutine assert_true
end program test_statistics
