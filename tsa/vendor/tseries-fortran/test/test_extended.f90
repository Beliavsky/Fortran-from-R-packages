! Experimental modern Fortran translation of computational routines from
! the R package tseries 0.10-62. Original authors include Adrian Trapletti
! and Kurt Hornik; Blake LeBaron contributed the original BDS code.
! Licensed under GPL-2.0-only OR GPL-3.0-only. See LICENSE and NOTICE.

program test_extended
   use tseries, only : dp, test_result, portfolio_result, seed_random, quadratic_map, pp_test, kpss_test, &
      adf_test, terasvirta_test, white_test, po_test, portfolio_optimize, fft_surrogate, amplitude_surrogate
   implicit none
   real(dp), allocatable :: x(:), y(:), z(:), mat(:, :), returns(:, :)
   type(test_result) :: r
   type(portfolio_result) :: pf
   integer :: n,i

   n=160
   call seed_random(44)
   x=quadratic_map(0.3_dp,3.8_dp,n)
   allocate(y(n),z(n),mat(n,2),returns(n,3))
   call fft_surrogate(x,y)
   call amplitude_surrogate(x,z)
   if(abs(sum(y)-sum(x))>1.0e-7_dp*max(1.0_dp,abs(sum(x)))) error stop 'fft surrogate mean'

   r=adf_test(x,lags=1)
   if(r%status/=0) error stop 'adf'
   r=pp_test(x,use_t_statistic=.true.)
   if(r%status/=0) error stop 'pp'
   r=kpss_test(x,trend=.true.)
   if(r%status/=0) error stop 'kpss'
   r=terasvirta_test(x,lag=2)
   if(r%status/=0) error stop 'terasvirta'
   r=white_test(x,lag=2,qstar=2,q=5,seed=22)
   if(r%status/=0) error stop 'white'

   mat(:,1)=x
   mat(:,2)=0.5_dp*x+0.05_dp*y
   r=po_test(mat)
   if(r%status/=0) error stop 'po'

   do i=1,n
      returns(i,1)=0.001_dp+0.02_dp*(x(i)-0.5_dp)
      returns(i,2)=0.0005_dp+0.01_dp*(y(i)-0.5_dp)
      returns(i,3)=0.0008_dp+0.015_dp*(z(i)-0.5_dp)
   end do
   pf=portfolio_optimize(returns,target_return=sum(returns)/real(size(returns),dp),shorts=.true.)
   if(pf%status/=0) error stop 'portfolio'
   if(abs(sum(pf%weights)-1.0_dp)>1.0e-6_dp) error stop 'portfolio sum'
   print '(a)','Extended tests passed.'
end program test_extended
