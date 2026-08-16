program test_rng
   use discretedists
   implicit none
   integer,parameter::n=12000
   integer::x(n),fails
   real(dp)::m
   fails=0
   call rcompo(n,2.0_dp,1.0_dp,x);m=sum(real(x,dp))/real(n,dp)
   call chk(abs(m-2.0_dp)<0.08_dp,'COMPO RNG mean')
   call rdgeii(n,0.4_dp,1.0_dp,x);m=sum(real(x,dp))/real(n,dp)
   call chk(abs(m-(0.4_dp/0.6_dp))<0.06_dp,'DGEII RNG mean')
   call rggeo(n,0.4_dp,1.0_dp,x);m=sum(real(x,dp))/real(n,dp)
   call chk(abs(m-(0.4_dp/0.6_dp))<0.06_dp,'GGEO RNG mean')
   call rhyperpo2(n,2.0_dp,1.0_dp,x);m=sum(real(x,dp))/real(n,dp)
   call chk(abs(m-2.0_dp)<0.08_dp,'HYPERPO2 RNG mean')
   call chk(minval(x)>=0,'RNG support')
   if(fails/=0)error stop 1
   print *,'test_rng: PASS'
contains
   subroutine chk(ok,name)
      logical,intent(in)::ok;character(*),intent(in)::name
      if(.not.ok)then;print *,trim(name),' failed';fails=fails+1;end if
   end subroutine
end program
