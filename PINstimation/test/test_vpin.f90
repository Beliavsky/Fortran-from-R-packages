! SPDX-License-Identifier: GPL-3.0-or-later
program test_vpin
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use pinstimation
   implicit none
   real(dp) :: changes(12),volume(12),duration(12),vbs,manual
   real(dp),allocatable :: bv(:),sv(:),bd(:),v(:),iv(:)
   integer :: status

   changes=[-1.0_dp,0.5_dp,1.2_dp,-0.2_dp,0.0_dp,0.8_dp,-1.1_dp,0.4_dp,1.0_dp,-0.7_dp,0.2_dp,0.3_dp]
   volume=[100.0_dp,120.0_dp,80.0_dp,110.0_dp,90.0_dp,130.0_dp,70.0_dp,100.0_dp,120.0_dp,80.0_dp,90.0_dp,110.0_dp]
   duration=60.0_dp
   call build_volume_buckets(changes,volume,duration,2,3,bv,sv,bd,vbs,status)
   if(status/=0.or.size(bv)/=6) error stop 'bucket construction failed'
   if(abs(vbs-200.0_dp)>1.0e-12_dp) error stop 'wrong VBS'
   if(abs(sum(bv+sv)-sum(volume))>1.0e-9_dp) error stop 'volume not conserved'
   call compute_vpin_from_buckets(bv,sv,2,v,bucket_size=vbs,status=status)
   manual=(abs(bv(1)-sv(1))+abs(bv(2)-sv(2)))/(2.0_dp*vbs)
   if(.not.ieee_is_nan(v(1))) error stop 'first VPIN should be missing'
   if(abs(v(2)-manual)>1.0e-12_dp) error stop 'wrong rolling VPIN'
   call compute_ivpin_from_buckets(bv,sv,bd,2,iv,status)
   if(status/=0.or.any(.not.ieee_is_nan(iv(2:)).and.iv(2:)<0.0_dp)) error stop 'IVPIN failed'
   if(any(iv(2:)>1.0_dp)) error stop 'IVPIN outside range'
   print '(a)', 'test_vpin: PASS'
end program test_vpin
