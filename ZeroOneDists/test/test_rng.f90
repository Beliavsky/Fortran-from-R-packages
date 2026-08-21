program test_rng
   use zero_one_dists
   implicit none
   integer,parameter::n=30000
   real(dp),allocatable::x(:),muv(:),sigv(:),nuv(:)
   real(dp)::m
   integer::fails,i
   allocate(x(n),muv(n),sigv(n),nuv(n))
   fails=0
   call rber2(x,0.63_dp,5.0_dp,0.4_dp)
   if(any(x<0.0_dp .or. x>1.0_dp))then
      print '(a)','BER2 RNG range FAIL'
      fails=fails+1
   end if
   m=sum(x)/real(n,dp)
   if(abs(m-0.63_dp)>0.012_dp)then
      print '(a,f12.6)','BER2 RNG mean FAIL ',m
      fails=fails+1
   end if
   call ruhlg(x,1.2_dp)
   if(any(x<=0.0_dp .or. x>=1.0_dp))then
      print '(a)','UHLG RNG range FAIL'
      fails=fails+1
   end if
   call rumb(x,0.8_dp)
   muv=0.6_dp+0.4_dp*[(real(i-1,dp)/real(n-1,dp),i=1,n)]
   call rumb(x,muv)
   if(any(x<=0.0_dp .or. x>=1.0_dp))then
      print '(a)','UMB vector-parameter RNG range FAIL'
      fails=fails+1
   end if
   if(any(x<=0.0_dp .or. x>=1.0_dp))then
      print '(a)','UMB RNG range FAIL'
      fails=fails+1
   end if
   call ruphn(x,1.4_dp,0.7_dp)
   if(any(x<=0.0_dp .or. x>=1.0_dp))then
      print '(a)','UPHN RNG range FAIL'
      fails=fails+1
   end if
   if(fails/=0)error stop 1
   print '(a)','test_rng: PASS'
end program test_rng
