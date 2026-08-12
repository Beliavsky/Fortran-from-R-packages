program school_matching
   use matchingmarkets
   implicit none
   integer :: sp(3,5),cp(5,3),slots(3)
   type(assignment_result_t) :: da
   integer :: s
   sp(:,1)=[1,2,3];sp(:,2)=[1,3,2];sp(:,3)=[2,1,3]
   sp(:,4)=[2,3,1];sp(:,5)=[3,2,1]
   cp(:,1)=[2,1,3,4,5];cp(:,2)=[4,3,1,2,5];cp(:,3)=[5,2,1,3,4]
   slots=[2,2,1]
   da=iaa(sp,cp,slots,'deferred')
   print '(a)', 'student -> college'
   do s=1,size(da%assignment)
      print '(i0,a,i0)',s,' -> ',da%assignment(s)
   end do
end program school_matching
