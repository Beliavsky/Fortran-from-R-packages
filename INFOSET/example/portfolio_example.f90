! SPDX-License-Identifier: GPL-2.0-or-later
program portfolio_example
  use infoset, only : dp, ptf_construction, portfolio_result
  use infoset, only : summary_ptf, portfolio_summary
  implicit none
  integer, parameter :: observations = 180, assets = 4
  real(dp) :: prices(observations,assets), r
  type(portfolio_result) :: portfolio
  type(portfolio_summary) :: statistics
  integer :: i, j
  prices(1,:) = [100.0_dp, 75.0_dp, 130.0_dp, 90.0_dp]
  do i = 2, observations
    do j = 1, assets
      r = 0.00015_dp*real(j,dp) + 0.004_dp*sin(0.08_dp*real(i*j,dp))
      prices(i,j) = prices(i-1,j)*exp(r)
    end do
  end do
  call ptf_construction(prices, 80, 20, 'M', portfolio)
  call summary_ptf(portfolio%oos_returns, statistics)
  write(*,'(a,*(f10.6,1x))') 'first-window weights: ', portfolio%weights(:,1)
  write(*,'(a,f10.6)') 'mean out-of-sample portfolio return: ', statistics%mean
end program portfolio_example
