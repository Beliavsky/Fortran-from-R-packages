program test_galeshapley
   use matchingr
   implicit none
   real(dp) :: um(2,3), uw(3,2)
   type(marriage_result_t) :: r
   type(college_result_t) :: c
   integer :: slots(2)
   um = reshape([0.0_dp,1.0_dp, 1.0_dp,0.0_dp, 0.0_dp,1.0_dp],[2,3])
   uw = reshape([0.0_dp,2.0_dp,1.0_dp, 1.0_dp,0.0_dp,2.0_dp],[3,2])
   r=marriage_market(um,uw)
   if(any(r%engagements /= [2,3])) error stop "marriage engagements"
   if(any(r%proposals /= [0,1,2])) error stop "marriage proposals"
   if(.not.gale_shapley_stable(um,uw,reshape(r%proposals,[3,1]),reshape(r%engagements,[2,1]))) then
      error stop "marriage stability"
   end if
   slots=[2,2]
   c=college_admissions(um,uw,slots,.true.)
   if(any(c%matched_students /= [2,1,2])) error stop "college students"
   if(c%matched_colleges(1,1)/=2 .or. c%matched_colleges(2,1)/=3 .or. c%matched_colleges(2,2)/=1) then
      error stop "college slots"
   end if
   print *, "test_galeshapley: PASS"
end program
