program test_estimation
   use discrete_weibull
   implicit none
   integer(i64) :: x1(20),x3(20)
   type(dweibull_fit_result) :: fp,fml,f3p,f3ml
   integer :: fails

   fails = 0
   x1 = [1_i64,1_i64,1_i64,1_i64,1_i64,2_i64,2_i64,2_i64,3_i64,3_i64, &
         3_i64,4_i64,4_i64,5_i64,5_i64,6_i64,7_i64,8_i64,9_i64,12_i64]
   fp = estdweibull(x1,"P")
   if (fp%status /= 0) fails=fails+1
   if (fp%pars(1) <= 0.0_dp .or. fp%pars(1) >= 1.0_dp .or. fp%pars(2) <= 0.0_dp) &
      fails=fails+1

   fml = estdweibull(x1,"ML")
   if (fml%status /= 0) fails=fails+1
   if (fml%objective > loglikedw([0.5_dp,1.0_dp],x1)+1.0e-8_dp) fails=fails+1

   x3 = [0_i64,0_i64,0_i64,0_i64,0_i64,1_i64,1_i64,1_i64,1_i64,1_i64, &
         1_i64,1_i64,2_i64,2_i64,2_i64,2_i64,3_i64,3_i64,4_i64,6_i64]
   f3p = estdweibull3(x3,"P")
   if (f3p%status /= 0) fails=fails+1
   f3ml = estdweibull3(x3,"ML")
   if (f3ml%status /= 0) fails=fails+1
   if (f3ml%objective > loglikedw3([1.0_dp,0.0_dp],x3)+1.0e-8_dp) fails=fails+1

   if (fails /= 0) then
      print *, "test_estimation: FAIL", fails
      print *, fp%pars, fml%pars, f3p%pars, f3ml%pars
      error stop 1
   end if
   print *, "test_estimation: PASS"
end program test_estimation
