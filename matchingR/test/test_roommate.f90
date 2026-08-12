program test_roommate
   use matchingr
   implicit none
   integer :: p1(3,4), p2(3,4), p3(2,3)
   type(roommate_result_t) :: r
   p1=reshape([2,3,4, 1,3,4, 1,2,4, 1,2,3],[3,4])
   r=stable_roommates_preferences(p1)
   if(.not.r%stable_exists) error stop "roommate p1 no match"
   if(any(r%matching /= [2,1,4,3])) error stop "roommate p1"
   if(.not.roommate_stable(p1,r%matching,.false.)) error stop "roommate p1 unstable"
   p2=reshape([4,3,2, 1,3,4, 2,1,4, 1,2,3],[3,4])
   r=stable_roommates_preferences(p2)
   if(.not.r%stable_exists) error stop "roommate p2 no match"
   if(any(r%matching /= [4,3,2,1])) error stop "roommate p2"
   p3=reshape([2,3, 1,3, 1,2],[2,3])
   r=stable_roommates_preferences(p3)
   if(.not.r%stable_exists) error stop "roommate odd no match"
   if(any(r%matching /= [2,1,0])) error stop "roommate odd"
   print *, "test_roommate: PASS"
end program
