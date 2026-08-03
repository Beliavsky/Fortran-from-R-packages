! SPDX-License-Identifier: GPL-2.0-or-later
program left_risk_example
  use infoset, only : dp, lr_cp, left_risk_result
  implicit none
  integer, parameter :: observations = 180, assets = 3
  real(dp) :: prices(observations,assets), r
  type(left_risk_result) :: result
  integer :: i, j
  prices(1,:) = [100.0_dp, 80.0_dp, 120.0_dp]
  do i = 2, observations
    do j = 1, assets
      r = 0.0003_dp + 0.003_dp*sin(0.11_dp*real(i*j,dp))
      if (mod(i+5*j,47) == 0) r = r - 0.04_dp
      prices(i,j) = prices(i-1,j)*exp(r)
    end do
  end do
  call lr_cp(prices, 80, 20, result)
  write(*,'(a,*(f10.6,1x))') 'first change points: ', result%first_change_point
  write(*,'(a,*(f10.6,1x))') 'first-window Left Risk: ', result%values(:,1)
end program left_risk_example
