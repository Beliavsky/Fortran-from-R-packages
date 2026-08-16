program test_bootstrap
   use powerlaw
   implicit none
   type(powerlaw_dist) :: m
   type(bootstrap_result) :: b,bp
   real(dp) :: x(18), xmins(4)
   integer :: i,fails
   fails=0
   do i=1,18
      x(i)=0.5_dp+0.25_dp*real(i,dp)
   end do
   xmins=[0.75_dp,1.0_dp,1.25_dp,1.5_dp]
   m=conexp(x)
   b=bootstrap(m,2,seed=123,xmax=6.0_dp,xmins=xmins)
   if(size(b%xmin)/=2 .or. b%successful<1) fails=fails+1
   bp=bootstrap_p(m,2,seed=123,xmax=6.0_dp,xmins=xmins)
   if(size(bp%xmin)/=2 .or. bp%p<0.0_dp .or. bp%p>1.0_dp) fails=fails+1
   if(fails/=0) then
      print *,"test_bootstrap: FAIL",fails
      error stop 1
   end if
   print *,"test_bootstrap: PASS"
end program
