program test_startvalues
   use gkwdist
   implicit none
   real(dp),allocatable :: x(:),par(:)
   real(dp) :: obj
   logical :: ok
   integer :: fails
   fails=0
   call seed_rng(777)
   x=rbeta_(80,2.5_dp,3.0_dp)
   par=gkwgetstartvalues(x,'beta',4,obj,ok)
   if(size(par)/=2 .or. any(par<=0.0_dp) .or. obj/=obj) then
      print '(a)','beta starting-value failure'; fails=fails+1
   end if
   par=gkwgetstartvalues(x,'kw',4,obj,ok)
   if(size(par)/=2 .or. any(par<=0.0_dp) .or. obj/=obj) then
      print '(a)','kw starting-value failure'; fails=fails+1
   end if
   if(fails==0) then
      print '(a)','test_startvalues: PASS'
   else
      print '(a,i0)','test_startvalues: FAIL ',fails
      error stop 1
   end if
end program test_startvalues
