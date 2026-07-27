! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_monte_carlo
   use optionpricing, only : dp, greeks_result, conditional_result, &
      asian_call_naive_mc, asian_call_ncv_lr_mc, asian_call_best_mc, &
      conditional_estimates_z, asian_call_app_lord
   implicit none
   type(greeks_result) :: naive1,naive2,ncv,best
   type(conditional_result) :: fixed
   real(dp) :: z(4,5),target
   integer :: i,j

   naive1=asian_call_naive_mc(3000,1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp,123)
   naive2=asian_call_naive_mc(3000,1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp,123)
   ncv=asian_call_ncv_lr_mc(3000,1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp,123)
   best=asian_call_best_mc(2400,1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp, &
      300,60,1.0e-13_dp,123)
   if(any([naive1%status,ncv%status,best%status]/=0)) error stop 1
   if(maxval(abs(naive1%estimate-naive2%estimate))>0.0_dp) error stop 1
   target=asian_call_app_lord(1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp,.true.)
   if(abs(best%estimate(1)-target)>2.0e-4_dp) error stop 1
   if(abs(best%estimate(2)-0.59385_dp)>2.0e-3_dp) error stop 1
   if(best%error95(1)>=ncv%error95(1)) error stop 1
   if(ncv%error95(1)>=naive1%error95(1)) error stop 1

   do j=1,5
      do i=1,4
         z(i,j)=sin(real(3*i+7*j,dp))
      end do
   end do
   fixed=conditional_estimates_z(z,1.0_dp,4,100.0_dp,0.05_dp,0.2_dp,100.0_dp, &
      'std',maxiter=100,tol=1.0e-14_dp)
   if(fixed%status/=0) error stop 1
   ! Transcendental-library implementations and Newton stopping can differ by
   ! several ulps across compilers. These limits still require about 10 to 11
   ! significant decimal digits for the fixed conditional-Monte-Carlo values.
   call assert_close(fixed%y(1),1.3997081977210567e-4_dp,1.0e-14_dp,1.0e-10_dp)
   call assert_close(fixed%y(2),-2.2543424968300073e-5_dp,1.0e-14_dp,1.0e-10_dp)
   call assert_close(fixed%controls(1),-2.0255613633713971e-4_dp,1.0e-14_dp,1.0e-10_dp)
   print '(a)', 'test_monte_carlo: PASS'
contains
   subroutine assert_close(x,y,atol,rtol)
      real(dp), intent(in) :: x,y,atol,rtol
      real(dp) :: err,limit
      err=abs(x-y)
      limit=max(atol,rtol*max(abs(x),abs(y)))
      if(err>limit) then
         print '(a,4(es24.16,1x))','mismatch: ',x,y,err,limit
         error stop 1
      end if
   end subroutine assert_close
end program test_monte_carlo
