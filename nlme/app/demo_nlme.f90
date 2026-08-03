! SPDX-License-Identifier: GPL-2.0-or-later
program demo_nlme
  use nlme
  implicit none
  integer,parameter :: ng=5,m=6,n=ng*m
  real(dp)::y(n),x(n,2),z(n,1),tt
  integer::group(n),g,i,k
  type(lme_result)::fit
  type(pd_spec)::random

  k=0
  do g=1,ng
    do i=1,m
      k=k+1;tt=real(i-1,dp)/5.0_dp
      x(k,:)=[1.0_dp,tt];z(k,1)=1.0_dp;group(k)=g
      y(k)=1.0_dp+2.0_dp*tt+0.3_dp*sin(real(g,dp))+0.05_dp*cos(real(k,dp))
    end do
  end do
  random%kind=PD_DIAG;random%dim=1;allocate(random%par(1));random%par=log(0.2_dp)
  call fit_lme(y,x,z,group,fit,random=random)
  write(*,'(a)')'nlme-fortran demo'
  write(*,'(a,2f11.5)')'fixed effects: ',fit%beta
  write(*,'(a,f11.5)')'random variance: ',fit%random_covariance(1,1)
  write(*,'(a,f11.5)')'residual sigma: ',fit%sigma
  write(*,'(a,f11.5)')'log likelihood: ',fit%log_likelihood
end program demo_nlme
