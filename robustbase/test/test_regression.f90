! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
program test_regression
   use robustbase
   use test_support
   implicit none
   integer,parameter::n=100
   real(dp)::x(n,2),y(n),u,pv
   type(robust_regression_result)::lts,mm,logit,pois
   type(robust_nls_result)::nls
   real(dp)::xn(60,1),yn(60),start(2)
   integer::i
   call seed_rng(24680)
   do i=1,n
      x(i,1)=1.0_dp;x(i,2)=-2.0_dp+4.0_dp*real(i-1,dp)/real(n-1,dp)
      y(i)=1.0_dp+2.0_dp*x(i,2)+0.08_dp*sin(real(i,dp))
   end do
   y(n-14:n)=y(n-14:n)+25.0_dp
   call lts_regression(x,y,lts,alpha=0.75_dp,n_starts=250)
   call mm_regression(x,y,mm,n_starts=200)
   call assert_true(abs(lts%coefficients(2)-2.0_dp)<0.15_dp,'LTS slope')
   call assert_true(abs(mm%coefficients(2)-2.0_dp)<0.15_dp,'MM slope')
   do i=1,n
      pv=1.0_dp/(1.0_dp+exp(-(-0.4_dp+1.3_dp*x(i,2))));call random_number(u);y(i)=merge(1.0_dp,0.0_dp,u<pv)
   end do
   y(n-4:n)=1.0_dp-y(n-4:n)
   call robust_glm_fit(x,y,'binomial',logit)
   call assert_true(logit%coefficients(2)>0.4_dp .and. logit%converged,'robust binomial')
   do i=1,n
      pv=exp(0.3_dp+0.35_dp*x(i,2));y(i)=real(poisson_draw(pv),dp)
   end do
   y(n-3:n)=y(n-3:n)+30.0_dp
   call robust_glm_fit(x,y,'poisson',pois)
   call assert_true(pois%coefficients(2)>0.05_dp .and. pois%converged,'robust poisson')
   do i=1,60
      xn(i,1)=2.0_dp*real(i-1,dp)/59.0_dp
      yn(i)=1.5_dp*exp(0.7_dp*xn(i,1))+0.03_dp*sin(real(i,dp))
   end do
   yn(56:60)=yn(56:60)+10.0_dp;start=[1.0_dp,0.3_dp]
   call robust_nls_fit(exp_model,xn,yn,start,nls,psi='huber',tuning=1.345_dp,max_iter=100)
   call assert_true(abs(nls%parameters(1)-1.5_dp)<0.35_dp .and. abs(nls%parameters(2)-0.7_dp)<0.25_dp,'robust nonlinear fit')
   write(*,'(a)')'LTS, MM, robust GLM, and nonlinear regression tests passed.'
contains
   integer function poisson_draw(lambda) result(k)
      real(dp),intent(in)::lambda
      real(dp)::l,prod,r
      l=exp(-lambda);prod=1.0_dp;k=0
      do
         k=k+1;call random_number(r);prod=prod*r
         if(prod<=l)exit
      end do
      k=k-1
   end function

end program test_regression
