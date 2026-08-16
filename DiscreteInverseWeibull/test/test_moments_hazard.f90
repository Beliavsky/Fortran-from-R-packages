program test_moments_hazard
 use discrete_inverse_weibull
 implicit none
 type(diw_moments)::e
 real(dp)::h,a
 integer::f
 f=0;e=ediweibull(0.5_dp,2.5_dp)
 if(abs(e%ex-1.7281293704708605_dp)>2e-13_dp)f=f+1
 if(abs(e%ex2-4.837443769558661_dp)>3e-13_dp)f=f+1
 h=hrdiweibull(1.0_dp,0.5_dp,2.5_dp)
 if(abs(h-0.5_dp)>1e-14_dp)f=f+1
 a=ahrdiweibull(1.0_dp,0.5_dp,2.5_dp)
 if(abs(a+log(0.5_dp))>1e-14_dp)f=f+1
 e=ediweibull(0.75_dp,1.25_dp)
 if(.not.(e%ex<huge(1.0_dp).and.e%ex2>1e300_dp))f=f+1
 if(f/=0)then;print *,'test_moments_hazard: FAIL',f;error stop 1;end if
 print *,'test_moments_hazard: PASS'
end program
