program test_preference_indexing
   use matchingr
   implicit none
   integer :: p1(2,3), p0(2,3), r1(3,2), r0(3,2)
   type(marriage_result_t) :: a,b
   p1=reshape([1,2, 2,1, 1,2],[2,3])
   p0=p1-1
   r1=reshape([1,2,3, 2,3,1],[3,2])
   r0=r1-1
   if(.not.check_preferences(p1)) error stop "1-based preferences"
   if(.not.check_preferences(p0)) error stop "0-based preferences"
   a=marriage_market_preferences(p1,r1)
   b=marriage_market_preferences(p0,r0)
   if(any(a%proposals/=b%proposals)) error stop "indexing mismatch"
   print *, "test_preference_indexing: PASS"
end program
