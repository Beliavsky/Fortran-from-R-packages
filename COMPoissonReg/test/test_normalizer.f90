program test_normalizer
   use compoissonreg
   implicit none
   type(cmp_control_t) :: ctrl
   real(dp) :: z1,z2
   integer :: t,fails
   fails=0
   ctrl=cmp_control_t()
   z1=ncmp(2.0_dp,1.0_dp,.true.,ctrl)
   if(abs(z1-2.0_dp)>5.0e-7_dp)then;print *,'FAIL logZ Poisson',z1;fails=fails+1;end if
   z2=cmp_logz_approx(100.0_dp,1.0_dp)
   if(abs(z2-100.0_dp)>1.0e-10_dp)then;print *,'FAIL approx Poisson',z2;fails=fails+1;end if
   t=tcmp(0.5_dp,1.0_dp,ctrl)
   if(t<5 .or. t>100)then;print *,'FAIL truncation',t;fails=fails+1;end if
   if(fails==0)then;print *,'test_normalizer: PASS';else;error stop 1;end if
end program test_normalizer
