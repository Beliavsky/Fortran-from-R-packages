program test_reference
   use discrete_laplace
   implicit none
   integer :: x(11),st,fails
   real(dp) :: ml(2),mm(2)
   fails=0
   x=[-3,-2,-1,-1,0,0,0,1,1,2,4]
   ml=estdlaplace2(x,'ML',status=st)
   if(st/=0) fails=fails+1
   if(abs(ml(1)-0.54984579_dp)>2.0e-6_dp) fails=fails+1
   if(abs(ml(2)-0.39785595_dp)>2.0e-6_dp) fails=fails+1
   if(abs(dlaplacelike2(ml,x)-22.397250840457904_dp)>2.0e-9_dp) fails=fails+1
   mm=estdlaplace2(x,'M',status=st)
   if(st/=0) fails=fails+1
   if(abs(mm(1)-0.52211356_dp)>2.0e-6_dp) fails=fails+1
   if(abs(mm(2)-0.34833350_dp)>2.0e-6_dp) fails=fails+1
   if(fails/=0) then
      print *, 'test_reference: FAIL',fails,ml,mm
      error stop 1
   end if
   print *, 'test_reference: PASS'
end program
