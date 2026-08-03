! SPDX-License-Identifier: GPL-3.0-only
program test_nuisance
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use esback
  implicit none
  integer, parameter :: n = 80
  real(dp) :: x(n,2), y(n), u(n)
  real(dp), allocatable :: mu(:), sigma(:), cv(:), dens(:), cdf(:), b(:)
  integer :: i, status

  do i = 1, n
    x(i,:) = [1.0_dp, sin(0.13_dp*i)]
    y(i) = 0.4_dp + 0.7_dp*x(i,2) + &
      (0.8_dp + 0.15_dp*x(i,2))*sin(1.91_dp*i)
  end do

  call conditional_mean_sigma(y, x, mu, sigma, status)
  call assert_true(status == esback_ok, 'location-scale status')
  call assert_true(all(sigma > 0.0_dp), 'positive scale')

  call quantile_regression(x, y, 0.1_dp, b, status)
  call assert_true(status == esback_ok, 'quantile regression status')
  u = y - matmul(x,b)

  call density_quantile_function(y, x, u, 0.1_dp, sparsity_nid, &
    bandwidth_hall_sheather, dens, status)
  call assert_true(status == esback_ok, 'NID density status')
  call assert_true(all(dens >= 0.0_dp), 'nonnegative density')

  call cdf_at_quantile(y, x, matmul(x,b), cdf, status)
  call assert_true(status == esback_ok, 'CDF status')
  call assert_true(all(cdf >= 0.0_dp .and. cdf <= 1.0_dp), 'CDF range')

  call conditional_truncated_variance(u, x, sigma_scl_n, cv, status)
  call assert_true(status == esback_ok, 'normal truncated variance status')
  call assert_true(all(ieee_is_finite(cv)) .and. all(cv >= 0.0_dp), 'valid truncated variance')

  print '(a)', 'test_nuisance: PASS'
contains
  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    if (.not. condition) then
      print *, trim(label)
      error stop 1
    end if
  end subroutine assert_true
end program test_nuisance
