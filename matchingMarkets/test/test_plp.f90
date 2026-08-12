program test_plp
   use matchingmarkets
   implicit none
   real(dp)::v(4,4)
   type(plp_result_t)::r
   v=0.0_dp
   v(1,2)=5;v(2,1)=5;v(3,4)=4;v(4,3)=4
   v(1,3)=1;v(3,1)=1;v(1,4)=1;v(4,1)=1
   v(2,3)=1;v(3,2)=1;v(2,4)=1;v(4,2)=1
   r=plp(v)
   if(abs(r%objective-18.0_dp)>1e-8_dp) error stop 'PLP objective'
   if(size(r%pairs,2)/=2) error stop 'PLP pair count'
   print *, 'test_plp: PASS'
end program
