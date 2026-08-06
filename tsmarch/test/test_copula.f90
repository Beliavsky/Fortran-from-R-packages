program test_copula
  use ghyp_kinds, only : dp
  use tsd_math, only : normal_cdf
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tsmarch
  use test_support
  implicit none
  type(copula_spec) :: spec
  type(dcc_parameters) :: p
  type(dcc_filter_result) :: filtered
  real(dp), allocatable :: z(:, :), u(:, :), recovered(:, :)
  integer :: i, n

  n = 120
  allocate(z(n,2),u(n,2))
  do i=1,n
    z(i,1) = sin(0.17_dp*real(i,dp)) + 0.3_dp*cos(0.07_dp*real(i,dp))
    z(i,2) = 0.55_dp*z(i,1) + 0.7_dp*cos(0.11_dp*real(i,dp))
    u(i,1) = min(max(normal_cdf(z(i,1)),1.0e-8_dp),1.0_dp-1.0e-8_dp)
    u(i,2) = min(max(normal_cdf(z(i,2)),1.0e-8_dp),1.0_dp-1.0e-8_dp)
  end do
  recovered = copula_transform(u,'gaussian')
  call assert_true(maxval(abs(recovered-z)) < 5.0e-6_dp, 'Gaussian PIT inverse')

  spec%distribution='gaussian'
  spec%alpha_order=1
  spec%gamma_order=0
  spec%beta_order=1
  allocate(p%alpha(1),p%gamma(0),p%beta(1))
  p%alpha=0.05_dp
  p%beta=0.90_dp
  p%shape=8.0_dp
  filtered = copula_filter(u,spec,p)
  call assert_true(filtered%status==tsm_success,'Gaussian copula filter')
  call assert_true(ieee_is_finite(filtered%log_likelihood),'Gaussian copula finite likelihood')

  spec%distribution='student'
  p%shape=7.0_dp
  filtered = copula_filter(u,spec,p)
  call assert_true(filtered%status==tsm_success,'Student copula filter')
  call assert_true(ieee_is_finite(filtered%log_likelihood),'Student copula finite likelihood')

  call finish_test('test_copula')
end program test_copula
