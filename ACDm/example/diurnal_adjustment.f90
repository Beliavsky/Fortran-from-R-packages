! SPDX-License-Identifier: GPL-3.0-or-later
program diurnal_adjustment
  use acdm
  implicit none
  integer, parameter :: n = 240
  real(dp) :: clock_time(n), durations(n)
  type(diurnal_result) :: result
  integer :: i

  do i = 1, n
    clock_time(i) = 36000.0_dp + real(mod(i - 1, 120), dp) * 4800.0_dp / 119.0_dp
    durations(i) = 2.0_dp + 0.7_dp * cos(2.0_dp * pi * &
      (clock_time(i) - 36000.0_dp) / 4800.0_dp) + 0.05_dp * sin(real(i, dp))
  end do

  call diurnal_adjust(clock_time, durations, result, DIURNAL_FFF, &
                      fourier_order=4)
  if (result%status /= ACDM_SUCCESS) error stop 'diurnal adjustment failed'
  print '(a,f10.5)', 'Mean raw duration:      ', sum(durations) / real(n, dp)
  print '(a,f10.5)', 'Mean adjusted duration: ', sum(result%adjusted) / real(n, dp)
  print '(a,3f10.5)', 'First adjusted values:  ', result%adjusted(1:3)
end program diurnal_adjustment
