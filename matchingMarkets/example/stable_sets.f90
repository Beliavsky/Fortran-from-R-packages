program stable_sets
   use matchingmarkets
   implicit none
   integer :: sp(3,3),cp(3,3),slots(3)
   type(stable_set_t)::sets
   integer::k
   sp(:,1)=[1,2,3];sp(:,2)=[2,1,3];sp(:,3)=[1,2,3]
   cp(:,1)=[2,1,3];cp(:,2)=[1,2,3];cp(:,3)=[1,2,3]
   slots=1
   sets=hri(sp,cp,slots,100)
   print '(a,i0)','stable matchings: ',sets%count
   do k=1,sets%count
      print '(*(i0,1x))',sets%assignments(:,k)
   end do
end program stable_sets
