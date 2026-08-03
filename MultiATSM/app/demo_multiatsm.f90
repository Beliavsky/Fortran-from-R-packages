program demo_multiatsm
  use multiatsm, only : dp, affine_loadings, response_result, variance_decomposition_result, &
    affine_yield_loadings, impulse_responses, forecast_error_variance_decomposition
  implicit none
  real(dp) :: phi(2, 2), sigma(2, 2), identity(2, 2), d(2)
  integer :: maturities(4), info
  type(affine_loadings) :: loadings
  type(response_result) :: irf
  type(variance_decomposition_result) :: fevd

  phi = reshape([0.92_dp, 0.02_dp, 0.04_dp, 0.80_dp], [2, 2])
  sigma = reshape([0.0004_dp, 0.00005_dp, 0.00005_dp, 0.0002_dp], [2, 2])
  identity = 0.0_dp
  identity(1, 1) = 1.0_dp
  identity(2, 2) = 1.0_dp
  d = [1.0_dp, 0.2_dp]
  maturities = [1, 12, 36, 60]
  call affine_yield_loadings(maturities, phi, 0.001_dp, sigma, d, loadings, info)
  if (info /= 0) error stop 'affine loadings failed'
  call impulse_responses(phi, sigma, identity, loadings%b, 12, irf, info)
  if (info /= 0) error stop 'IRF failed'
  call forecast_error_variance_decomposition(phi, sigma, identity, loadings%b, 12, fevd, info)
  if (info /= 0) error stop 'FEVD failed'
  write(*, '(a,4f11.6)') 'Affine intercepts: ', loadings%a
  write(*, '(a,2f11.6)') '60-month yield impact responses: ', irf%yields(4, :, 1)
  write(*, '(a,2f11.6)') '60-month yield FEVD at 12 periods: ', fevd%yields(4, :, 12)
end program demo_multiatsm
