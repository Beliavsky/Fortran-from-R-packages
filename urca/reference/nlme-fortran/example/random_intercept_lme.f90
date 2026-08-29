! SPDX-License-Identifier: GPL-2.0-or-later
program random_intercept_lme
  use nlme
  implicit none
  integer, parameter :: ng=8,m=6,n=ng*m
  real(dp) :: y(n),x(n,2),z(n,1),time(n),b,tt
  integer :: group(n),g,i,k
  type(pd_spec) :: random
  type(lme_result) :: fit
  type(nlme_control) :: control

  k=0
  do g=1,ng
    b=0.4_dp*sin(real(g,dp))
    do i=1,m
      k=k+1
      tt=real(i-1,dp)/real(m-1,dp)
      x(k,:)=[1.0_dp,tt]
      z(k,1)=1.0_dp
      group(k)=g
      time(k)=real(i,dp)
      y(k)=2.0_dp+1.1_dp*tt+b+0.06_dp*cos(1.3_dp*real(k,dp))
    end do
  end do
  random%kind=PD_DIAG
  random%dim=1
  allocate(random%par(1))
  random%par=log(0.25_dp)
  control%max_iter=500
  call fit_lme(y,x,z,group,fit,random=random,time=time,control=control)
  write(*,'(a,2f12.6)')'fixed effects: ',fit%beta
  write(*,'(a,f12.6)')'random-intercept variance: ',fit%random_covariance(1,1)
  write(*,'(a,f12.6)')'residual sigma: ',fit%sigma
  write(*,'(a)')'BLUPs:'
  do g=1,size(fit%group_levels)
    write(*,'(i5,f12.6)')fit%group_levels(g),fit%random_effects(g,1)
  end do
end program random_intercept_lme
