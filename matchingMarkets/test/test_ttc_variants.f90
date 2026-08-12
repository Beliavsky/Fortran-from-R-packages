program test_ttc_variants
   use matchingmarkets
   implicit none
   integer :: sp(3,3),cp(3,3),slots(3),prio(3),kp(4,3)
   type(assignment_result_t)::a,b
   sp(:,1)=[1,2,3];sp(:,2)=[2,1,3];sp(:,3)=[3,1,2]
   cp(:,1)=[2,1,3];cp(:,2)=[1,2,3];cp(:,3)=[3,1,2]
   slots=1;prio=[1,2,3]
   a=ttc_school(sp,cp,slots,prio)
   if(count(a%assignment>0)/=3) error stop 'ttc_school matching size'
   kp(:,1)=[2,4,1,3];kp(:,2)=[1,4,2,3];kp(:,3)=[4,1,2,3]
   b=ttcc_kidney(kp,prio)
   if(any(b%assignment<0)) error stop 'ttcc invalid assignment'
   print *, 'test_ttc_variants: PASS'
end program
