! SPDX-License-Identifier: GPL-2.0-or-later
program nonlinear_gnls_example
  use nlme
  implicit none
  integer, parameter :: n=40
  real(dp) :: x(n,1),y(n),theta0(3)
  integer :: i
  type(nonlinear_result) :: fit
  type(nlme_control) :: control

  do i=1,n
    x(i,1)=real(i-1,dp)/15.0_dp
  end do
  y=2.4_dp*exp(-0.8_dp*x(:,1))+0.3_dp+0.01_dp*sin([(real(i,dp),i=1,n)])
  theta0=[1.8_dp,0.5_dp,0.1_dp]
  control%max_iter=200
  call fit_gnls(exponential_decay_model,y,x,theta0,fit,control=control)
  write(*,'(a,3f12.6)')'parameters: ',fit%parameters
  write(*,'(a,f12.6)')'sigma: ',fit%sigma
end program nonlinear_gnls_example
