program test_eadam
   use matchingmarkets
   implicit none
   integer :: sp(3,4),cp(4,3),slots(3)
   logical :: consent(4)
   type(eadam_result_t)::ea
   type(assignment_result_t)::da
   sp(:,1)=[1,2,3];sp(:,2)=[1,2,3];sp(:,3)=[2,1,3];sp(:,4)=[2,3,1]
   cp(:,1)=[2,1,3,4];cp(:,2)=[4,3,1,2];cp(:,3)=[1,2,3,4]
   slots=[1,1,1];consent=.false.
   ea=hri3_eadam(sp,cp,slots,'deferred',consent,100)
   da=iaa(sp,cp,slots,'deferred')
   if(any(ea%assignment/=da%assignment)) error stop 'EADAM no-consent differs from DA'
   if(ea%ea_iterations/=0) error stop 'EADAM iteration count'
   print *, 'test_eadam: PASS'
end program
