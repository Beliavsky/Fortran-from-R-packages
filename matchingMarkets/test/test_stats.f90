program test_stats
   use matchingmarkets
   implicit none
   integer,parameter::n=120
   real(dp)::x(n,3),y(n)
   type(khb_result_t)::r
   integer::i
   do i=1,n
      x(i,1)=1.0_dp
      x(i,2)=real(i-n/2,dp)/20.0_dp
      x(i,3)=0.5_dp*x(i,2)+0.2_dp*sin(real(i,dp))
      y(i)=merge(1.0_dp,0.0_dp,0.7_dp*x(i,2)+0.5_dp*x(i,3)+0.8_dp*sin(1.7_dp*real(i,dp))>0.0_dp)
   end do
   r=khb(x,y,3)
   if(size(r%p_value)/=2) error stop 'KHB size'
   if(any(r%p_value<0.0_dp).or.any(r%p_value>1.0_dp)) error stop 'KHB pvalues'
   print *, 'test_stats: PASS'
end program
