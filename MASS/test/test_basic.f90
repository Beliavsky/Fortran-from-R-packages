! SPDX-License-Identifier: GPL-3.0-only
program test_basic
  use mass
  use test_support
  implicit none
  real(dp) :: a(2,2), integral, sigma(2,2), mu(2), empirical_cov(2,2)
  real(dp), allocatable :: inverse(:,:), basis(:,:), samples(:,:), rv(:)
  real(dp), allocatable :: contrasts(:,:)
  integer :: status
  type(kde2d_result) :: density
  integer(kind=selected_int_kind(18)), allocatable :: num(:), den(:)

  a = reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2,2])
  inverse = ginv(a, status=status)
  call assert_true(status == mass_success, 'ginv status')
  call assert_close(inverse(1,1), 2.0_dp/3.0_dp, 1.0e-10_dp, 'ginv fixture')
  call assert_close(inverse(1,2), -1.0_dp/3.0_dp, 1.0e-10_dp, 'ginv off diagonal')

  call null_space(reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp], [2,3]), &
    basis, status)
  call assert_true(status == mass_success .and. size(basis,2) == 1, 'null-space dimension')
  call assert_close(abs(basis(3,1)), 1.0_dp, 1.0e-10_dp, 'null-space basis')

  mu = [1.0_dp, -2.0_dp]
  sigma = reshape([1.0_dp,0.4_dp,0.4_dp,2.0_dp], [2,2])
  call mvrnorm(30, mu, sigma, samples, status, empirical=.true., seed=12345)
  call assert_true(status == mass_success, 'mvrnorm status')
  call assert_close(sum(samples(:,1))/30.0_dp, mu(1), 1.0e-10_dp, 'empirical mean 1')
  call assert_close(sum(samples(:,2))/30.0_dp, mu(2), 1.0e-10_dp, 'empirical mean 2')
  empirical_cov = matmul(transpose(samples-spread(mu,1,30)), &
    samples-spread(mu,1,30))/29.0_dp
  call assert_close(empirical_cov(1,1), sigma(1,1), 1.0e-8_dp, 'empirical covariance')

  call kde2d(samples(:,1), samples(:,2), density, n_grid=[12,10])
  call assert_true(density%status == mass_success, 'kde2d status')
  call assert_true(all(density%density >= 0.0_dp), 'kde2d nonnegative')

  integral = adaptive_area(square, 0.0_dp, 1.0_dp, limit=20, eps=1.0e-10_dp, status=status)
  call assert_close(integral, 1.0_dp/3.0_dp, 1.0e-8_dp, 'adaptive integration')

  call rational_approximation([0.5_dp, 1.0_dp/3.0_dp], num, den)
  call assert_true(num(1) == 1 .and. den(1) == 2, 'rational one half')
  call assert_true(num(2) == 1 .and. den(2) == 3, 'rational one third')
  rv = rational_values([0.5_dp, 1.0_dp/3.0_dp])
  call assert_close(rv(2), 1.0_dp/3.0_dp, 1.0e-12_dp, 'rational value')

  contrasts = successive_difference_contrasts(4)
  call assert_true(size(contrasts,1) == 4 .and. size(contrasts,2) == 3, 'contrast shape')
  call assert_true(all(abs(sum(contrasts,dim=1)) < 1.0e-12_dp), 'contrast columns sum zero')
  write(*,'(a)') 'test_basic: PASS'
contains
  function square(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value
    value = x*x
  end function square
end program test_basic
