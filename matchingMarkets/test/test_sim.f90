program test_sim
   use matchingmarkets
   implicit none
   type(one_sided_sim_t)::a
   type(two_sided_sim_t)::b
   a=stabsim(3,4,2,123_i8)
   if(size(a%market)/=24) error stop 'stabsim size'
   b=stabsim2(2,5,2,[2,2],456_i8)
   if(size(b%market)/=20) error stop 'stabsim2 size'
   if(count(b%matched==1)>8) error stop 'stabsim2 capacity'
   print *, 'test_sim: PASS'
end program
