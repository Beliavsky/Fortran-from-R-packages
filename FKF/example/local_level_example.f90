! SPDX-License-Identifier: GPL-2.0-or-later
program local_level_example
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use fkf_module
  implicit none

  integer, parameter :: n = 20
  type(fkf_model) :: model
  type(fkf_result) :: filtered
  type(fks_result) :: smoothed
  real(dp) :: y(1, n), level
  integer :: i

  model%a0 = [0.0_dp]
  model%p0 = reshape([10.0_dp], [1, 1])
  model%dt = reshape([0.0_dp], [1, 1])
  model%ct = reshape([0.0_dp], [1, 1])
  model%tt = reshape([1.0_dp], [1, 1, 1])
  model%zt = reshape([1.0_dp], [1, 1, 1])
  model%hht = reshape([0.08_dp], [1, 1, 1])
  model%ggt = reshape([0.25_dp], [1, 1, 1])

  level = 1.0_dp
  do i = 1, n
    level = level + 0.12_dp * sin(0.4_dp * real(i, dp))
    y(1, i) = level + 0.18_dp * cos(0.7_dp * real(i, dp))
  end do
  y(1, 7) = ieee_value(0.0_dp, ieee_quiet_nan)
  y(1, 14) = ieee_value(0.0_dp, ieee_quiet_nan)

  call kalman_filter(model, y, filtered, .true.)
  call kalman_smooth(model, y, filtered, smoothed)
  if (filtered%status /= fkf_success .or. smoothed%status /= fkf_success) error stop 1

  write(*, '(a,f12.6)') 'log likelihood: ', filtered%log_likelihood
  write(*, '(a)') ' t       observation       filtered        smoothed'
  do i = 1, n
    write(*, '(i2,3(2x,f14.6))') i, y(1, i), filtered%att(1, i), smoothed%ahatt(1, i)
  end do
end program local_level_example
