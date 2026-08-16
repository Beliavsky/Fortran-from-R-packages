program test_fit
   use powerlaw
   implicit none
   type(powerlaw_dist) :: m
   type(estimate_pars_result) :: e
   real(dp) :: x(100), target
   integer :: fails
   fails=0
   m=conpl([1.0_dp,2.0_dp,4.0_dp,8.0_dp])
   call m%set_xmin(1.0_dp)
   e=estimate_pars(m)
   target=1.0_dp+4.0_dp/log(64.0_dp)
   if(abs(e%pars(1)-target)>1.0e-10_dp) fails=fails+1

   m=conexp([1.0_dp,2.0_dp,3.0_dp,4.0_dp])
   call m%set_xmin(0.0_dp)
   e=estimate_pars(m)
   if(abs(e%pars(1)-0.4_dp)>1.0e-12_dp) fails=fails+1

   m=disexp([1.0_dp,2.0_dp,3.0_dp,4.0_dp])
   call m%set_xmin(1.0_dp)
   e=estimate_pars(m)
   target=log(1.0_dp+1.0_dp/1.5_dp)
   if(abs(e%pars(1)-target)>1.0e-12_dp) fails=fails+1

   x=2.0_dp
   m=dispois(x)
   call m%set_xmin(0.0_dp)
   e=estimate_pars(m)
   if(abs(e%pars(1)-2.0_dp)>3.0e-5_dp) fails=fails+1

   m=conlnorm([exp(-1.0_dp),1.0_dp,exp(1.0_dp)])
   call m%set_xmin(1.0e-12_dp)
   e=estimate_pars(m)
   if(abs(e%pars(1))>5.0e-4_dp) fails=fails+1
   if(abs(e%pars(2)-sqrt(2.0_dp/3.0_dp))>5.0e-4_dp) fails=fails+1

   if(fails/=0) then
      print *,"test_fit: FAIL",fails
      error stop 1
   end if
   print *,"test_fit: PASS"
end program
