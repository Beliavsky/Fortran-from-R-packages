! Experimental modern Fortran translation of computational routines from
! the R package tseries 0.10-62. Original authors include Adrian Trapletti
! and Kurt Hornik; Blake LeBaron contributed the original BDS code.
! Licensed under GPL-2.0-only OR GPL-3.0-only. See LICENSE and NOTICE.

program test_smoke
   use tseries, only : dp, test_result, arma_result, garch_result, bds_result, &
      quadratic_map, jarque_bera_test, runs_test, arma_fit, garch_fit, bds_test, &
      stationary_bootstrap, seed_random
   implicit none
   real(dp), allocatable :: x(:),boot(:)
   integer :: binary(8)
   type(test_result) :: test
   type(arma_result) :: ar
   type(garch_result) :: ga
   type(bds_result) :: bds

   x=quadratic_map(0.2_dp,4.0_dp,128)
   if(size(x)/=128) error stop 'quadratic_map size failed'
   if(abs(x(1)-0.2_dp)>1.0e-14_dp) error stop 'quadratic_map initial value failed'

   test=jarque_bera_test(x)
   if(test%status/=0 .or. test%statistic<0.0_dp) error stop 'Jarque-Bera failed'

   binary=[0,0,1,1,0,1,0,1]
   test=runs_test(binary)
   if(test%status/=0 .or. test%p_value<0.0_dp .or. test%p_value>1.0_dp) error stop 'runs test failed'

   ar=arma_fit(x,1,0,max_iterations=300)
   if(.not.allocated(ar%coefficients)) error stop 'ARMA allocation failed'
   if(ar%css<=0.0_dp) error stop 'ARMA CSS failed'

   ga=garch_fit(x-sum(x)/real(size(x),dp),1,1,max_iterations=300)
   if(.not.allocated(ga%coefficients)) error stop 'GARCH allocation failed'
   if(any(ga%coefficients<=0.0_dp)) error stop 'GARCH positivity failed'
   if(sum(ga%coefficients(2:))>=1.0_dp) error stop 'GARCH stationarity failed'

   bds=bds_test(x,max_embedding=3,eps=[0.1_dp])
   if(bds%status/=0) error stop 'BDS failed'

   allocate(boot(size(x))); call seed_random(777)
   call stationary_bootstrap(x,5.0_dp,boot)
   if(any(boot<minval(x)) .or. any(boot>maxval(x))) error stop 'bootstrap failed'

   print '(a)','All tests passed.'
end program test_smoke
