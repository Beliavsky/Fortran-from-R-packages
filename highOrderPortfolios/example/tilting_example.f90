! SPDX-License-Identifier: GPL-3.0-only
program tilting_example
   use highorderportfolios
   implicit none
   integer, parameter :: t=180,n=4
   real(dp) :: x(t,n),w0(n),base_moments(4),d(4),kappa
   type(sample_moments) :: statistics
   type(portfolio_result) :: solution
   integer :: i

   do i=1,t
      x(i,1)=0.0014_dp+0.011_dp*sin(0.09_dp*i)+0.003_dp*sin(0.023_dp*i)**2
      x(i,2)=0.0008_dp+0.008_dp*cos(0.12_dp*i)
      x(i,3)=0.0003_dp+0.013_dp*sin(0.065_dp*i+0.7_dp)
      x(i,4)=0.0010_dp+0.006_dp*cos(0.04_dp*i)+0.002_dp*sin(0.15_dp*i)**3
   end do

   call estimate_sample_moments(x,statistics,adjust_magnitude=.true.)
   w0=1.0_dp/real(n,dp)
   base_moments=eval_portfolio_moments(w0,statistics)
   d=abs(base_moments)
   kappa=0.30_dp*sqrt(dot_product(w0,matmul(statistics%covariance,w0)))
   call design_mvsktilting_portfolio_via_sample_moments(d,statistics,solution, &
      w_init=w0,w0=w0,w0_moments=base_moments,kappa=kappa,method='L-MVSKT')

   print '(a,*(f10.6,1x))','weights: ',solution%w
   print '(a,f10.6)','tilting delta: ',solution%delta
   print '(a,*(f10.6,1x))','relative improvements: ',solution%improvement
end program tilting_example
