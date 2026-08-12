program test_stable_sets
   use matchingmarkets
   implicit none
   integer :: sp(3,3),cp(3,3),slots(3),rp(4,4)
   type(stable_set_t)::hs,rs
   sp(:,1)=[1,2,3];sp(:,2)=[2,1,3];sp(:,3)=[1,2,3]
   cp(:,1)=[2,1,3];cp(:,2)=[1,2,3];cp(:,3)=[1,2,3]
   slots=1
   hs=hri_all(sp,cp,slots,100)
   if(hs%count<1) error stop 'hri_all no solution'
   rp(:,1)=[2,3,4,0];rp(:,2)=[1,3,4,0];rp(:,3)=[4,1,2,0];rp(:,4)=[3,1,2,0]
   rs=sri_all(rp,100)
   if(rs%count<1) error stop 'sri_all no solution'
   print *, 'test_stable_sets: PASS',hs%count,rs%count
end program
