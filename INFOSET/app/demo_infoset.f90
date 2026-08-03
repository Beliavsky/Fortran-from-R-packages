! SPDX-License-Identifier: GPL-2.0-or-later
program demo_infoset
  use infoset, only : dp, g_ret, infoset_estimate, information_set_result
  implicit none
  integer, parameter :: observations = 240
  real(dp) :: prices(observations), r
  real(dp), allocatable :: gross(:)
  type(information_set_result) :: fit
  integer :: i
  prices(1) = 100.0_dp
  do i = 2, observations
    r = 0.0004_dp + 0.005_dp*sin(0.17_dp*real(i,dp))
    if (mod(i,43) == 0) r = r - 0.055_dp
    prices(i) = prices(i-1)*exp(r)
  end do
  gross = g_ret(prices)
  call infoset_estimate(gross, fit)
  write(*,'(a,i0)') 'number of detected change points: ', fit%n_change_points
  if (fit%n_change_points > 0) then
    write(*,'(a,*(f10.6,1x))') 'change points: ', &
      fit%change_points(1:fit%n_change_points)
  end if
end program demo_infoset
