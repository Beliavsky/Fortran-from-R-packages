! SPDX-License-Identifier: GPL-3.0-only
program sample_mvsk_example
   use highorderportfolios
   implicit none
   integer, parameter :: t=180,n=4
   real(dp) :: x(t,n),lambda(4)
   type(sample_moments) :: moments
   type(portfolio_result) :: solution
   integer :: i

   do i=1,t
      x(i,1)=0.0012_dp+0.010_dp*sin(0.13_dp*i)+0.002_dp*sin(0.037_dp*i)**2
      x(i,2)=0.0007_dp+0.007_dp*cos(0.11_dp*i)
      x(i,3)=0.0004_dp+0.012_dp*sin(0.07_dp*i+0.8_dp)
      x(i,4)=0.0009_dp+0.006_dp*cos(0.05_dp*i)+0.001_dp*sin(0.17_dp*i)**3
   end do

   call estimate_sample_moments(x,moments)
   lambda=[1.0_dp,5.0_dp,18.333333333333333_dp,55.0_dp]
   call design_mvsk_portfolio_via_sample_moments(lambda,moments,solution, &
      method='Q-MVSK',maxiter=200)

   print '(a,*(f10.6,1x))','weights: ',solution%w
   print '(a,*(es13.5,1x))','moments: ',solution%moments
   print '(a,l1)','converged: ',solution%converged
end program sample_mvsk_example
