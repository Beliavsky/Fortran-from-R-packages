program affine_pricing_example
  use multiatsm, only : dp, affine_loadings, affine_yield_loadings, forecast_result, forecast_yields
  implicit none
  integer :: maturities(6), info, j
  real(dp) :: k1q(2, 2), sigma(2, 2), short_rate_loadings(2), k0(2), last_state(2)
  type(affine_loadings) :: loadings
  type(forecast_result) :: forecast

  maturities = [1, 3, 6, 12, 24, 60]
  k1q = reshape([0.96_dp, 0.0_dp, 0.03_dp, 0.85_dp], [2, 2])
  sigma = reshape([0.0004_dp, 0.00005_dp, 0.00005_dp, 0.0002_dp], [2, 2])
  short_rate_loadings = [1.0_dp, 0.25_dp]
  call affine_yield_loadings(maturities, k1q, 0.001_dp, sigma, short_rate_loadings, loadings, info)
  if (info /= 0) error stop 'affine_yield_loadings failed'
  write(*, '(a)') 'Maturity, intercept, first loading'
  do j = 1, size(maturities)
    write(*, '(i4,2f14.7)') maturities(j), loadings%a(j), loadings%b(j, 1)
  end do
  k0 = [0.001_dp, 0.0_dp]
  last_state = [0.025_dp, -0.01_dp]
  call forecast_yields(k0, k1q, loadings%a, loadings%b, last_state, 3, forecast, info)
  if (info /= 0) error stop 'forecast_yields failed'
  write(*, '(a,6f10.5)') 'One-step yields: ', forecast%yields(:, 1)
end program affine_pricing_example
