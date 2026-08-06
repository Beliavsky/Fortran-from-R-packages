! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
program test_statistics_denoise
  use waveslim
  use waveslim_test_support
  implicit none
  real(dp) :: clean(256), noisy(256)
  real(dp), allocatable :: reconstructed(:), spectrum(:), tapers(:,:)
  real(dp), allocatable :: gains(:), cascade(:)
  type(wavelet_transform) :: wt, shrunk
  type(interval_vector) :: variance_result, covariance_result
  integer :: i
  real(dp) :: mse_noisy, mse_denoised

  do i = 1, size(clean)
    clean(i) = sin(0.05_dp*real(i,dp)) + 0.5_dp*cos(0.16_dp*real(i,dp))
    noisy(i) = clean(i) + 0.25_dp*sin(1.91_dp*real(i,dp)) + &
      0.18_dp*cos(2.37_dp*real(i,dp))
  end do
  wt = dwt(noisy, 'la8', 5)
  shrunk = universal_thresh(wt, 4, .false.)
  reconstructed = idwt(shrunk)
  mse_noisy = sum((noisy-clean)**2)/real(size(clean),dp)
  mse_denoised = sum((reconstructed-clean)**2)/real(size(clean),dp)
  call assert_true(mse_denoised < mse_noisy, 'wavelet denoising improves MSE')

  variance_result = wave_variance(wt)
  covariance_result = wave_covariance(wt, wt)
  call assert_true(size(variance_result%estimate) == 6, &
    'wave variance scale count')
  call assert_close_scalar(maxval(abs(variance_result%estimate- &
    covariance_result%estimate)), 0.0_dp, 1.0e-10_dp, &
    'variance equals self covariance')

  spectrum = periodogram(clean)
  call assert_true(size(spectrum) == 129, 'periodogram size')
  call assert_true(all(spectrum >= 0.0_dp), 'periodogram nonnegative')

  tapers = sine_taper(64, 4)
  call assert_close_scalar(sum(tapers(:,1)**2), 1.0_dp, 1.0e-12_dp, &
    'sine taper normalization')

  cascade = wavelet_filter('haar', 'HL')
  call assert_true(size(cascade) == 4, 'cascade-filter size')
  gains = squared_gain('haar', 'HL', 128)
  call assert_true(size(gains) == 65, 'squared-gain size')
  call assert_true(minval(gains) >= -1.0e-12_dp, &
    'squared gain nonnegative')

  write(*,'(a)') 'test_statistics_denoise: PASS'
end program test_statistics_denoise
