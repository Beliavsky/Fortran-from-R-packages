! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
program test_transforms_1d
  use waveslim
  use waveslim_test_support
  implicit none
  real(dp) :: x(128), sum_parts(128)
  real(dp), allocatable :: y(:), y1(:), y2(:)
  type(wavelet_transform) :: wt
  type(mra_result) :: mr
  type(complex_wavelet_transform) :: hwt, dt
  type(real_vector), allocatable :: coh(:), phase(:)
  type(wavelet_filter_type) :: filter
  type(status_type) :: status
  integer :: i, j

  do i = 1, size(x)
    x(i) = sin(0.11_dp*real(i,dp)) + 0.35_dp*cos(0.37_dp*real(i,dp))
  end do

  filter = wave_filter('la8', status)
  call assert_true(status%ok(), 'LA8 filter lookup')
  call assert_true(filter%length() == 8, 'LA8 filter length')
  call assert_close_scalar(sum(filter%lpf), sqrt2, 1.0e-10_dp, &
    'LA8 low-pass normalization')

  wt = dwt(x, 'la8', 4)
  call assert_true(wt%status%ok(), 'DWT status')
  y = idwt(wt)
  call assert_close_array(y, x, 2.0e-10_dp, 'DWT inverse')

  wt = modwt(x, 'la8', 4)
  call assert_true(wt%status%ok(), 'MODWT status')
  y = imodwt(wt)
  call assert_close_array(y, x, 2.0e-10_dp, 'MODWT inverse')

  mr = mra(x, 'la8', 4, 'modwt')
  sum_parts = mr%smooth
  do j = 1, size(mr%detail)
    sum_parts = sum_parts + mr%detail(j)%values
  end do
  call assert_close_array(sum_parts, x, 5.0e-10_dp, 'MRA additive identity')

  hwt = dwt_hilbert(x, 'k3l3', 3)
  call assert_true(hwt%status%ok(), 'Hilbert DWT status')
  y1 = idwt_hilbert(hwt, 1)
  y2 = idwt_hilbert(hwt, 2)
  call assert_close_array(y1, x, 2.0e-8_dp, 'Hilbert tree 1 inverse')
  call assert_close_array(y2, x, 2.0e-8_dp, 'Hilbert tree 2 inverse')

  hwt = modwt_hilbert(x, 'k3l3', 3)
  coh = modhwt_coherence(hwt, hwt, 5)
  phase = modhwt_phase(hwt, hwt, 5)
  do j = 1, size(coh)
    call assert_true(minval(coh(j)%values) > 0.999999_dp, &
      'self coherence is one')
    call assert_true(maxval(abs(phase(j)%values)) < 1.0e-12_dp, &
      'self phase is zero')
  end do

  dt = dualtree(x, 3)
  call assert_true(dt%status%ok(), 'dual-tree status')
  y = idualtree(dt)
  call assert_close_array(y, x, 2.0e-8_dp, 'dual-tree inverse')

  write(*,'(a)') 'test_transforms_1d: PASS'
end program test_transforms_1d
