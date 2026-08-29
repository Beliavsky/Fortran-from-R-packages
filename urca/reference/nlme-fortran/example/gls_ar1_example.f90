! SPDX-License-Identifier: GPL-2.0-or-later
program gls_ar1_example
  use nlme
  implicit none
  integer, parameter :: n=50
  real(dp) :: y(n),x(n,2),time(n),e(n),innovation
  integer :: i
  type(correlation_spec) :: corr
  type(gls_result) :: fit
  type(nlme_control) :: control

  x(:,1)=1.0_dp
  time=[(real(i,dp),i=1,n)]
  x(:,2)=(time-25.5_dp)/25.0_dp
  e(1)=0.05_dp
  do i=2,n
    innovation=0.08_dp*sin(1.9_dp*real(i,dp))
    e(i)=0.7_dp*e(i-1)+innovation
  end do
  y=1.2_dp+0.6_dp*x(:,2)+e
  corr%kind=COR_AR1
  allocate(corr%par(1))
  corr%par=0.2_dp
  control%max_iter=400
  call fit_gls(y,x,fit,correlation=corr,time=time,control=control)
  write(*,'(a,2f12.6)')'fixed effects: ',fit%beta
  write(*,'(a,f12.6)')'AR(1) correlation: ',fit%correlation_parameters(1)
  write(*,'(a,f12.6)')'residual sigma: ',fit%sigma
end program gls_ar1_example
