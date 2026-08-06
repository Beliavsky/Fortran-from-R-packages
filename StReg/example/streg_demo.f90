! SPDX-License-Identifier: GPL-2.0-only
program streg_demo
   use streg, only : dp, streg_fit, streg_options, stlm
   implicit none
   integer, parameter :: n=80
   real(dp) :: x(n,1),y(n)
   integer :: i
   type(streg_fit) :: fit
   type(streg_options) :: options

   do i=1,n
      x(i,1)=-2.0_dp+4.0_dp*real(i-1,dp)/real(n-1,dp)
      y(i)=1.2_dp+2.3_dp*x(i,1)+0.25_dp*sin(1.7_dp*real(i,dp))
   end do
   options%max_iter=250
   fit=stlm(y,x,v=8.0_dp,options=options)

   write(*,'(a,l1)')'converged: ',fit%converged
   write(*,'(a,f12.6)')'log likelihood: ',fit%log_likelihood
   write(*,'(a,2f12.6)')'intercept and slope: ',fit%beta(1,:)
   write(*,'(a,f12.6)')'R squared: ',fit%r_squared(1)
   write(*,'(a,2f12.6)')'AD statistic and p value: ',fit%ad_test(1,:)
end program streg_demo
