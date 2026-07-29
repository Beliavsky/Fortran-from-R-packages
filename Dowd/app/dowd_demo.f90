! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program dowd_demo
  use dowd, only: dp, normal_var, normal_es, student_t_var, student_t_es, &
       historical_var, historical_es, black_scholes_call_price, &
       variance_covariance_var, variance_covariance_es
  implicit none

  real(dp), parameter :: cl = 0.99_dp
  real(dp), parameter :: hp = 10.0_dp
  real(dp) :: profit_loss(10)
  real(dp) :: covariance(2,2), mu(2), positions(2)

  profit_loss = [0.8_dp,-1.1_dp,0.3_dp,-2.5_dp,1.4_dp,-0.7_dp,0.2_dp,-3.2_dp,1.1_dp,-0.4_dp]
  covariance = reshape([0.0004_dp,0.0001_dp,0.0001_dp,0.0009_dp],[2,2])
  mu = [0.0005_dp,0.0003_dp]
  positions = [1.0_dp,0.5_dp]

  write(*,'(a,f12.6)') "Normal VaR:       ",normal_var(0.001_dp,0.02_dp,cl,hp)
  write(*,'(a,f12.6)') "Normal ES:        ",normal_es(0.001_dp,0.02_dp,cl,hp)
  write(*,'(a,f12.6)') "Student-t VaR:    ",student_t_var(0.001_dp,0.02_dp,6.0_dp,cl,hp)
  write(*,'(a,f12.6)') "Student-t ES:     ",student_t_es(0.001_dp,0.02_dp,6.0_dp,cl,hp)
  write(*,'(a,f12.6)') "Historical VaR:   ",historical_var(profit_loss,0.95_dp)
  write(*,'(a,f12.6)') "Historical ES:    ",historical_es(profit_loss,0.95_dp)
  write(*,'(a,f12.6)') "Portfolio VaR:    ",variance_covariance_var(covariance,mu,positions,0.99_dp,1.0_dp)
  write(*,'(a,f12.6)') "Portfolio ES:     ",variance_covariance_es(covariance,mu,positions,0.99_dp,1.0_dp)
  write(*,'(a,f12.6)') "BS call price:    ",black_scholes_call_price(100.0_dp,105.0_dp,0.03_dp,0.20_dp,0.5_dp)
end program dowd_demo
