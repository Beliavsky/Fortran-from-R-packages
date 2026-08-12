program test_mechanisms
   use matchingmarkets
   implicit none
   integer :: sp(2,4), cp(4,2), slots(2), prio(4), hp(4,4), houses(4)
   type(assignment_result_t) :: da,im,r,t
   integer, allocatable :: bp(:,:)
   sp(:,1)=[1,2]; sp(:,2)=[1,2]; sp(:,3)=[2,1]; sp(:,4)=[2,1]
   cp(:,1)=[2,1,3,4]; cp(:,2)=[4,3,1,2]
   slots=[1,1]
   da=iaa(sp,cp,slots,'deferred')
   im=iaa(sp,cp,slots,'immediate')
   if(any(da%assignment /= [0,1,0,2])) error stop 'DA result'
   if(any(im%assignment /= [0,1,0,2])) error stop 'IAA result'
   bp=stability_check(da%assignment,cp,sp,slots)
   if(size(bp,2)/=0) error stop 'DA stability'
   prio=[1,2,3,4]
   r=rsd(sp,[1,1],prio)
   if(r%assignment(1)/=1 .or. r%assignment(2)/=2) error stop 'RSD'
   hp(:,1)=[1,2,3,4];hp(:,2)=[2,1,3,4];hp(:,3)=[3,1,2,4];hp(:,4)=[4,1,2,3]
   houses=[1,2,3,4]
   t=ttc_tenants(hp,houses,prio)
   if(any(t%assignment<1)) error stop 'TTC tenant'
   print *, 'test_mechanisms: PASS'
end program
