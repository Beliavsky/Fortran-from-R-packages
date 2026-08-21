! SPDX-License-Identifier: GPL-3.0-only
! Derived from the GPL-3 R package poilog by Vidar Grotan and Steinar Engen.
program test_rng
   use poilog, only : dp,poilog_seed,rpoilog,rbipoilog
   implicit none
   integer, allocatable :: x(:),xy(:,:)
   call poilog_seed(12345)
   x=rpoilog(100,0.0_dp,1.0_dp,keep0=.true.)
   if(size(x)/=100 .or. any(x<0)) error stop 1
   x=rpoilog(50,0.0_dp,1.0_dp,cond_s=.true.)
   if(size(x)/=50 .or. any(x<=0)) error stop 1
   xy=rbipoilog(40,0.0_dp,0.0_dp,1.0_dp,1.0_dp,0.5_dp,cond_s=.true.)
   if(size(xy,1)/=40 .or. size(xy,2)/=2 .or. any(xy(:,1)+xy(:,2)<=0)) error stop 1
   print *, 'test_rng: PASS'
end program test_rng
