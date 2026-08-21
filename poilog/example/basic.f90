! SPDX-License-Identifier: GPL-3.0-only
! Derived from the GPL-3 R package poilog by Vidar Grotan and Steinar Engen.
program basic
   use poilog, only : dp,dpoilog,dbipoilog,poilog_seed,rpoilog
   implicit none
   integer, allocatable :: x(:)
   integer :: n
   do n=0,10
      write(*,'(i3,1x,es16.8)') n,dpoilog(n,0.0_dp,1.0_dp)
   end do
   write(*,'(a,es16.8)') 'P(N1=2,N2=3) = ',dbipoilog(2,3,0.0_dp,0.0_dp,1.0_dp,1.0_dp,0.5_dp)
   call poilog_seed(20260819)
   x=rpoilog(10,0.0_dp,1.0_dp,cond_s=.true.)
   write(*,'(a,*(i0,1x))') 'draws: ',x
end program basic
