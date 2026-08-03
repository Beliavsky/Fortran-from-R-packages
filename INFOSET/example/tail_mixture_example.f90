! SPDX-License-Identifier: GPL-2.0-or-later
program tail_mixture_example
  use infoset, only : dp, tail_mixture, tail_mixture_result
  implicit none
  real(dp) :: gross_returns(300)
  type(tail_mixture_result) :: fit
  integer :: i
  do i = 1, 120
    gross_returns(i) = exp(-0.07_dp + 0.02_dp*sin(0.31_dp*real(i,dp)))
  end do
  do i = 121, 300
    gross_returns(i) = exp(0.025_dp + 0.018_dp*sin(0.23_dp*real(i,dp)))
  end do
  call tail_mixture(gross_returns, 0.0_dp, 1, fit)
  write(*,'(a,f10.6)') 'change point: ', fit%change_point
  write(*,'(a,2f10.5)') 'component means: ', fit%left_mean, fit%right_mean
  write(*,'(a,2f10.5)') 'component standard deviations: ', fit%left_sd, fit%right_sd
  write(*,'(a,f10.5)') 'left prior probability: ', fit%left_probability
end program tail_mixture_example
