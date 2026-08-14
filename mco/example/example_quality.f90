! SPDX-License-Identifier: GPL-2.0-only
program example_quality
   use mco, only : dp, dominated_hypervolume, generational_distance, epsilon_indicator
   implicit none
   real(dp) :: front(2,3), approximate(2,3)
   front=reshape([1.0_dp,3.0_dp,2.0_dp,2.0_dp,3.0_dp,1.0_dp],[2,3])
   approximate=front+0.1_dp
   print '(a,f8.4)', 'Hypervolume: ',dominated_hypervolume(front,[4.0_dp,4.0_dp])
   print '(a,f8.4)', 'Generational distance: ',generational_distance(approximate,front)
   print '(a,f8.4)', 'Additive epsilon indicator: ',epsilon_indicator(approximate,front)
end program
