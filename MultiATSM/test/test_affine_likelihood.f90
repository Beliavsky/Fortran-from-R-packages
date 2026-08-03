program test_affine_likelihood
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use multiatsm_kinds, only : dp
  use multiatsm_types, only : affine_loadings, atsm_likelihood_result
  use multiatsm_affine, only : affine_yield_loadings, pricing_factor_loadings, build_yield_intercepts
  use multiatsm_likelihood, only : atsm_log_likelihood, yield_error_variance
  implicit none
  integer, parameter :: nt = 80
  integer :: maturities(3), info, t
  real(dp) :: k1q(1, 1), sigma(1, 1), d(1), latent_b(2, 1), wpca(1, 2)
  real(dp) :: states(1, nt), spanned(1, nt), yields(2, nt), we(1, 2)
  real(dp) :: k0(1), k1(1, 1), sigma_states(1, 1), a(2), b(2, 1)
  real(dp) :: wpca_full(2, 2), measurement(1)
  real(dp), allocatable :: observed_b(:, :), rotation(:, :), observed_a(:), variances(:), covariance(:, :)
  type(affine_loadings) :: loadings
  type(atsm_likelihood_result) :: likelihood

  maturities = [1, 2, 3]
  k1q(1, 1) = 0.8_dp
  sigma(1, 1) = 0.04_dp
  d = 1.0_dp
  call affine_yield_loadings(maturities, k1q, 0.002_dp, sigma, d, loadings, info)
  call check(info == 0, 'affine loading status')
  call check(abs(loadings%b(1, 1) - 1.0_dp) < 1.0e-12_dp, 'B1 loading')
  call check(abs(loadings%b(2, 1) - 0.9_dp) < 1.0e-12_dp, 'B2 loading')
  call check(abs(loadings%b(3, 1) - 2.44_dp / 3.0_dp) < 1.0e-12_dp, 'B3 loading')
  call check(abs(loadings%a(2) + 0.008_dp) < 1.0e-12_dp, 'A2 loading')

  latent_b = reshape([1.0_dp, 0.5_dp], [2, 1])
  wpca = reshape([0.6_dp, 0.8_dp], [1, 2])
  call pricing_factor_loadings(latent_b, wpca, observed_b, rotation, info)
  call check(info == 0, 'pricing rotation status')
  call check(abs(rotation(1, 1) - 1.0_dp) < 1.0e-12_dp, 'pricing rotation')
  call build_yield_intercepts([0.01_dp, 0.02_dp], observed_b, wpca, observed_a, info)
  call check(info == 0, 'intercept construction')
  call check(abs(dot_product(wpca(1, :), observed_a)) < 1.0e-12_dp, 'PCA intercept restriction')

  k0 = 0.01_dp
  k1(1, 1) = 0.7_dp
  states(1, 1) = 0.02_dp
  do t = 2, nt
    states(1, t) = k0(1) + k1(1, 1) * states(1, t - 1) + 0.004_dp * sin(0.4_dp * real(t, dp))
  end do
  spanned = states
  a = [0.01_dp, 0.015_dp]
  b(:, 1) = [0.8_dp, 1.2_dp]
  yields = spread(a, 2, nt) + matmul(b, spanned)
  yields(2, :) = yields(2, :) + [(0.001_dp * cos(0.3_dp * real(t, dp)), t = 1, nt)]
  we = reshape([-0.8_dp, 0.6_dp], [1, 2])
  sigma_states(1, 1) = 0.004_dp**2 / 2.0_dp
  call atsm_log_likelihood(yields, states, spanned, wpca, we, a, b, k0, k1, sigma_states, 1, &
    likelihood, info)
  call check(info == 0, 'likelihood status')
  call check(ieee_is_finite(likelihood%negative_log_likelihood), 'finite likelihood')
  call check(likelihood%yield_error_variance(1) > 0.0_dp, 'measurement variance')

  wpca_full = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
  measurement = 0.25_dp
  call yield_error_variance(wpca_full, 1, measurement, 1, variances, covariance, info)
  call check(info == 0, 'yield covariance status')
  call check(abs(variances(1)) < 1.0e-12_dp .and. abs(variances(2) - 0.25_dp) < 1.0e-12_dp, &
    'yield covariance mapping')
  print '(a)', 'test_affine_likelihood: PASS'
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // message
      error stop 1
    end if
  end subroutine check
end program test_affine_likelihood
