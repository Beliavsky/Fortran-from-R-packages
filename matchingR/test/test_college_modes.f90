program test_college_modes
   use matchingr
   implicit none
   real(dp) :: su(3,4), cu(4,3)
   integer :: slots(3)
   integer :: prop(4,1), eng(3,1)
   type(college_result_t) :: a,b
   su=reshape([ &
      0.9_dp,0.5_dp,0.2_dp, &
      0.8_dp,0.7_dp,0.1_dp, &
      0.3_dp,0.9_dp,0.6_dp, &
      0.4_dp,0.2_dp,0.95_dp],[3,4])
   cu=reshape([ &
      0.9_dp,0.7_dp,0.5_dp,0.3_dp, &
      0.2_dp,0.95_dp,0.6_dp,0.4_dp, &
      0.8_dp,0.1_dp,0.7_dp,0.9_dp],[4,3])
   slots=[1,1,1]
   a=college_admissions(su,cu,slots,.true.)
   b=college_admissions(su,cu,slots,.false.)
   prop(:,1)=a%matched_students
   eng(:,1)=a%matched_colleges(:,1)
   if(.not.gale_shapley_stable(su,cu,prop,eng)) error stop "student-optimal unstable"
   prop(:,1)=b%matched_students
   eng(:,1)=b%matched_colleges(:,1)
   if(.not.gale_shapley_stable(su,cu,prop,eng)) error stop "college-optimal unstable"
   if(count(a%matched_students>0)/=3 .or. count(b%matched_students>0)/=3) error stop "wrong cardinality"
   print *, "test_college_modes: PASS"
end program test_college_modes
