program test_helpers
   use matchingmarkets
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none
   integer,allocatable::p(:,:),c(:,:)
   real(dp)::s(2,5,2)
   real(dp),allocatable::t(:,:)
   integer::i
   p=pair_combinations([1,3,5])
   if(size(p,2)/=3)error stop 'pairs'
   c=coalition_partitions(2,2,.false.)
   if(size(c,2)/=6)error stop 'coalitions'
   do i=1,5
      s(:,i,1)=[real(i,dp),real(i*i,dp)]
      s(:,i,2)=[real(i,dp)+0.1_dp,real(i*i+i,dp)-0.1_dp]
   end do
   t=consensus_mc(s)
   if(any(ieee_is_nan(t)))error stop 'consensus nan'
   print *, 'test_helpers: PASS'
end program
