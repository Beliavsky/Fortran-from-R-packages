program test_couples
   use matchingmarkets
   implicit none
   integer :: sp(2,4),cp(4,2),slots(2),members(2,1),choice(2,2,1)
   type(assignment_result_t)::r
   sp(:,1)=[1,2];sp(:,2)=[2,1];sp(:,3)=[1,2];sp(:,4)=[2,1]
   cp(:,1)=[1,3,2,4];cp(:,2)=[2,4,1,3]
   slots=[2,2];members(:,1)=[3,4]
   choice(1,:,1)=[1,2];choice(2,:,1)=[2,1]
   r=hri2_couples_exact(sp,cp,slots,members,choice,100000)
   if(size(r%assignment)/=4) error stop 'hri2 size'
   if(count(r%assignment>0)<2) error stop 'hri2 no useful matching'
   print *, 'test_couples: PASS'
end program
