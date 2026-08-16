program test_reference
   use discrete_weibull
   implicit none
   integer(i64) :: x1(20),x3(20)
   type(dweibull_fit_result) :: f1,f3
   integer :: fails

   fails = 0
   x1 = [1_i64,1_i64,1_i64,1_i64,1_i64,2_i64,2_i64,2_i64,3_i64,3_i64, &
         3_i64,4_i64,4_i64,5_i64,5_i64,6_i64,7_i64,8_i64,9_i64,12_i64]
   f1 = estdweibull(x1,"ML")
   if (abs(f1%pars(1)-0.7786074436731096_dp) > 2.0e-6_dp) fails=fails+1
   if (abs(f1%pars(2)-1.0843725801480233_dp) > 2.0e-6_dp) fails=fails+1
   if (abs(f1%objective-44.91058679693818_dp) > 2.0e-7_dp) fails=fails+1

   x3 = [0_i64,0_i64,0_i64,0_i64,0_i64,1_i64,1_i64,1_i64,1_i64,1_i64, &
         1_i64,1_i64,2_i64,2_i64,2_i64,2_i64,3_i64,3_i64,4_i64,6_i64]
   f3 = estdweibull3(x3,"ML")
   if (abs(f3%pars(1)-0.3355878526045979_dp) > 2.0e-6_dp) fails=fails+1
   if (abs(f3%pars(2)-0.5880910565736708_dp) > 2.0e-6_dp) fails=fails+1
   if (abs(f3%objective-32.96591337480952_dp) > 2.0e-7_dp) fails=fails+1

   if (fails /= 0) then
      print *, "test_reference: FAIL", fails, f1%pars, f3%pars
      error stop 1
   end if
   print *, "test_reference: PASS"
end program test_reference
