! SPDX-License-Identifier: GPL-3.0-only
program skew_t_mvsk_example
   use highorderportfolios
   implicit none
   integer, parameter :: t=80,n=3
   real(dp) :: x(t,n),lambda(4)
   type(skew_t_parameters) :: parameters
   type(portfolio_result) :: solution
   integer :: i

   do i=1,t
      x(i,1)=0.018_dp*sin(0.11_dp*i)+0.004_dp*sin(0.029_dp*i)**2
      x(i,2)=0.014_dp*cos(0.08_dp*i)+0.002_dp*sin(0.043_dp*i)**3
      x(i,3)=0.016_dp*sin(0.06_dp*i+0.6_dp)+0.003_dp*cos(0.021_dp*i)**2
   end do

   call estimate_skew_t(x,parameters,nu_lb=9.0_dp,max_iter=30,pxem=.true.)
   if(parameters%status/=hop_success .and. parameters%status/=hop_not_converged) then
      print '(a)','fit failed: '//trim(parameters%message)
      stop 1
   end if
   lambda=[1.0_dp,4.0_dp,10.0_dp,20.0_dp]
   call design_mvsk_portfolio_via_skew_t(lambda,parameters,solution, &
      method='PGD',maxiter=150,initial_eta=10.0_dp)

   print '(a,f10.4)','estimated nu: ',parameters%nu
   print '(a,*(f10.6,1x))','weights: ',solution%w
   print '(a,*(es13.5,1x))','moments: ',solution%moments
end program skew_t_mvsk_example
