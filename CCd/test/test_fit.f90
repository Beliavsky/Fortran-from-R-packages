program test_fit
   use ccd, only : dp, cc_fit_result, cc_reg_result, loc0_test_result, &
      cc_mle, cc_mle0, cc_reg, loc0_test
   implicit none
   real(dp) :: y(9), yr(9), x(9,1)
   type(cc_fit_result) :: f0, f1
   type(cc_reg_result) :: fr
   type(loc0_test_result) :: lt
   integer :: fail, i
   fail = 0
   y = real([-4,-2,-1,0,0,1,1,3,5], dp)
   f0 = cc_mle0(y, 1.0e-10_dp)
   f1 = cc_mle(y)
   call check_close(f0%lambda, 1.3114468924283338_dp, 2.0e-6_dp, fail)
   call check_close(f0%loglik, -22.229588534495807_dp, 2.0e-8_dp, fail)
   call check_close(f1%mu, 0.20469188_dp, 2.0e-5_dp, fail)
   call check_close(f1%lambda, 1.2918989368_dp, 2.0e-5_dp, fail)
   call check_close(f1%loglik, -22.17581012147_dp, 2.0e-7_dp, fail)

   x(:,1) = [(-2.0_dp + 0.5_dp*real(i-1,dp), i=1,9)]
   yr = real([-2,-2,-1,0,0,1,1,3,3], dp)
   fr = cc_reg(yr, x)
   call check_close(fr%lambda, 0.4329481214_dp, 2.0e-4_dp, fail)
   call check_close(fr%beta(1), 0.33603203_dp, 2.0e-4_dp, fail)
   call check_close(fr%beta(2), 1.31467036_dp, 2.0e-4_dp, fail)
   call check_close(fr%loglik, -8.2367562785_dp, 2.0e-5_dp, fail)

   lt = loc0_test(y)
   call check_close(lt%statistic, 0.10755682605_dp, 5.0e-5_dp, fail)
   if (lt%p_value <= 0.0_dp .or. lt%p_value >= 1.0_dp) fail = fail + 1

   if (fail == 0) then
      print *, 'test_fit: PASS'
   else
      print *, 'test_fit: FAIL', fail
      error stop 1
   end if
contains
   subroutine check_close(x, ref, eps, fail)
      real(dp), intent(in) :: x, ref, eps
      integer, intent(inout) :: fail
      if (abs(x-ref) > eps) then
         print *, 'mismatch:', x, ref, 'diff=', x-ref
         fail = fail + 1
      end if
   end subroutine check_close
end program test_fit
