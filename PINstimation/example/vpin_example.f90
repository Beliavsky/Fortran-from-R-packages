! SPDX-License-Identifier: GPL-3.0-or-later
program vpin_example
   use pinstimation
   implicit none
   real(dp) :: changes(12),volume(12),duration(12)
   type(vpin_result) :: estimate
   integer :: status

   changes=[-1.0_dp,0.5_dp,1.2_dp,-0.2_dp,0.0_dp,0.8_dp,-1.1_dp,0.4_dp,1.0_dp,-0.7_dp,0.2_dp,0.3_dp]
   volume=[100.0_dp,120.0_dp,80.0_dp,110.0_dp,90.0_dp,130.0_dp,70.0_dp,100.0_dp,120.0_dp,80.0_dp,90.0_dp,110.0_dp]
   duration=60.0_dp
   call compute_vpin(changes,volume,duration,2,3,2,estimate,improved=.true.,status=status)
   print '(a,f10.2)', 'volume bucket size: ',estimate%volume_bucket_size
   print '(a,*(f10.5,1x))', 'VPIN:  ',estimate%vpin
   print '(a,*(f10.5,1x))', 'IVPIN: ',estimate%ivpin
end program vpin_example
