! SPDX-License-Identifier: GPL-3.0-only
program var_t_example
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use imputefin
  implicit none
  integer,parameter::n=20
  real(dp)::y(n,2)
  type(var_t_result)::fit
  type(var_t_options)::opt
  integer::i
  y(1,:)=[0.0_dp,0.0_dp]
  do i=2,n
    y(i,1)=0.1_dp+0.6_dp*y(i-1,1)+0.1_dp*y(i-1,2)+0.03_dp*sin(real(i,dp))
    y(i,2)=-0.05_dp+0.2_dp*y(i-1,1)+0.4_dp*y(i-1,2)+0.02_dp*cos(real(i,dp))
  end do
  y(8,1)=ieee_value(0.0_dp,ieee_quiet_nan)
  opt%p=1
  call fit_var_t(y,fit,opt)
  write(*,'(a,f10.4)')'nu: ',fit%nu
  write(*,'(a,2f10.4)')'intercept: ',fit%phi0
  write(*,'(a)')'lag-1 matrix:'
  do i=1,2;write(*,'(2f12.6)')fit%phi(i,:,1);end do
end program var_t_example
