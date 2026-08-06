! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
program waveslim_demo
  use waveslim
  implicit none
  integer, parameter :: n = 256
  real(dp) :: clean(n), noisy(n)
  real(dp), allocatable :: reconstructed(:), denoised(:)
  type(wavelet_transform) :: wt, thresholded
  integer :: i
  real(dp) :: reconstruction_error, noisy_rmse, denoised_rmse

  do i = 1, n
    if (i <= 64) then
      clean(i) = 0.0_dp
    else if (i <= 128) then
      clean(i) = 2.0_dp
    else if (i <= 192) then
      clean(i) = -1.0_dp
    else
      clean(i) = 1.0_dp
    end if
    noisy(i) = clean(i) + 0.45_dp*(sin(1.731_dp*real(i,dp)) &
      + 0.7_dp*cos(2.437_dp*real(i,dp)))
  end do

  wt = dwt(noisy, 'la8', 5)
  if (.not. wt%status%ok()) error stop trim(wt%status%message)
  reconstructed = idwt(wt)
  thresholded = sure_thresh(wt, 4, .false.)
  denoised = idwt(thresholded)

  reconstruction_error = maxval(abs(reconstructed-noisy))
  noisy_rmse = sqrt(sum((noisy-clean)**2)/real(n,dp))
  denoised_rmse = sqrt(sum((denoised-clean)**2)/real(n,dp))

  write(*,'(a,es12.4)') 'DWT reconstruction maximum error: ', reconstruction_error
  write(*,'(a,es12.4)') 'Noisy RMSE:                    ', noisy_rmse
  write(*,'(a,es12.4)') 'SURE-thresholded RMSE:         ', denoised_rmse
end program waveslim_demo
