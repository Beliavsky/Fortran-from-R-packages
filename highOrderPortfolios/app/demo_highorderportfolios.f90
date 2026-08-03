! SPDX-License-Identifier: GPL-3.0-only
program demo_highorderportfolios
   use highorderportfolios
   implicit none
   integer, parameter :: t=120,n=3
   real(dp) :: x(t,n),lambda(4)
   type(sample_moments) :: statistics
   type(portfolio_result) :: result
   integer :: i

   do i=1,t
      x(i,1)=0.001_dp+0.010_dp*sin(0.10_dp*i)
      x(i,2)=0.0007_dp+0.007_dp*cos(0.08_dp*i)
      x(i,3)=0.0009_dp+0.009_dp*sin(0.06_dp*i+0.5_dp)
   end do
   call estimate_sample_moments(x,statistics)
   lambda=[1.0_dp,4.0_dp,10.0_dp,20.0_dp]
   call design_mvsk_portfolio_via_sample_moments(lambda,statistics,result,method='MM')
   print '(a,*(f9.5,1x))','MVSK weights: ',result%w
end program demo_highorderportfolios
