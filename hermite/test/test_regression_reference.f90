program test_regression_reference
   use hermite
   implicit none
   integer(i64) :: y(30)
   real(dp) :: x(30,1)
   type(hermite_glm_result) :: flog,fid
   integer :: fails

   fails=0
   y=[0_i64,1_i64,0_i64,2_i64,3_i64,0_i64,1_i64,4_i64,3_i64,6_i64, &
      0_i64,2_i64,1_i64,3_i64,5_i64,0_i64,1_i64,2_i64,3_i64,7_i64, &
      0_i64,4_i64,1_i64,3_i64,6_i64,2_i64,0_i64,1_i64,5_i64,3_i64]
   x(:,1)=1.0_dp

   ! Independent SciPy optimization of the exact m=3 likelihood gives
   ! mu=2.2999999860, d=2.1401904374, NLL=59.070751508180344.
   flog=fit_glm_hermite(y,x,link=HERMITE_LINK_LOG,m=3)
   if (flog%status>1) fails=fails+1
   if (abs(exp(flog%beta(1))-2.299999986012609_dp) > 2.0e-5_dp) fails=fails+1
   if (abs(flog%dispersion-2.1401904374149785_dp) > 3.0e-5_dp) fails=fails+1
   if (abs(flog%loglik+59.070751508180344_dp) > 3.0e-8_dp) fails=fails+1

   fid=fit_glm_hermite(y,x,link=HERMITE_LINK_IDENTITY,m=3, &
                       start_beta=[2.0_dp],start_d=2.0_dp)
   if (fid%status>1) fails=fails+1
   if (abs(fid%beta(1)-2.299999986012609_dp) > 2.0e-5_dp) fails=fails+1
   if (abs(fid%dispersion-2.1401904374149785_dp) > 3.0e-5_dp) fails=fails+1

   if (fails/=0) then
      print *,'test_regression_reference: FAIL',fails
      print *,'log:',flog%beta,flog%dispersion,flog%loglik
      print *,'id:',fid%beta,fid%dispersion,fid%loglik
      error stop 1
   end if
   print *,'test_regression_reference: PASS'
end program test_regression_reference
