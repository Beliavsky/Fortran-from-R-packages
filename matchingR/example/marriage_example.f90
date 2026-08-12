program marriage_example
   use matchingr
   implicit none
   real(dp) :: um(2,3), uw(3,2)
   type(marriage_result_t) :: res
   um=reshape([0.0_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,3])
   uw=reshape([0.0_dp,2.0_dp,1.0_dp,1.0_dp,0.0_dp,2.0_dp],[3,2])
   res=marriage_market(um,uw)
   print '(a,*(i0,1x))', 'proposer -> reviewer: ',res%proposals
   print '(a,*(i0,1x))', 'reviewer -> proposer: ',res%engagements
end program
