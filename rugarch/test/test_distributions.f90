! Part of the experimental modern Fortran translation of rugarch 1.5-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original rugarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-3.0-only

program test_distributions
   use rugarch
   implicit none
   real(dp), parameter :: probs(5)=[0.01_dp,0.10_dp,0.50_dp,0.90_dp,0.99_dp]
   real(dp), parameter :: tail_probs(7)=[1.0e-10_dp,1.0e-6_dp,0.01_dp,0.5_dp, &
      0.99_dp,1.0_dp-1.0e-6_dp,1.0_dp-1.0e-10_dp]
   real(dp), parameter :: t_shapes(4)=[2.1_dp,3.0_dp,7.0_dp,30.0_dp]
   real(dp), parameter :: ged_shapes(4)=[0.5_dp,1.0_dp,1.6_dp,4.0_dp]
   real(dp) :: q, p
   integer :: i, j

   do i=1,size(probs)
      q=qstd(probs(i),0.0_dp,1.0_dp,7.0_dp)
      p=pstd(q,0.0_dp,1.0_dp,7.0_dp)
      call assert_close(p,probs(i),2.0e-6_dp,'std round trip')
      q=qsged(probs(i),0.0_dp,1.0_dp,1.6_dp,1.3_dp)
      p=psged(q,0.0_dp,1.0_dp,1.6_dp,1.3_dp)
      call assert_close(p,probs(i),3.0e-6_dp,'sged round trip')
      q=qjsu(probs(i),0.0_dp,1.0_dp,0.5_dp,1.8_dp)
      p=pjsu(q,0.0_dp,1.0_dp,0.5_dp,1.8_dp)
      call assert_close(p,probs(i),2.0e-7_dp,'jsu round trip')
   end do
   do j=1,size(t_shapes)
      do i=1,size(tail_probs)
         q=qstd(tail_probs(i),0.0_dp,1.0_dp,t_shapes(j))
         p=pstd(q,0.0_dp,1.0_dp,t_shapes(j))
         call assert_close(p,tail_probs(i),2.0e-9_dp,'std accelerated inverse tails')
      end do
   end do
   do j=1,size(ged_shapes)
      do i=1,size(tail_probs)
         q=qged(tail_probs(i),0.0_dp,1.0_dp,ged_shapes(j))
         p=pged(q,0.0_dp,1.0_dp,ged_shapes(j))
         call assert_close(p,tail_probs(i),2.0e-9_dp,'GED accelerated inverse tails')
      end do
   end do
   print '(a)', 'distribution tests passed'
contains
   subroutine assert_close(a,b,tol,label)
      real(dp), intent(in) :: a,b,tol
      character(len=*), intent(in) :: label
      if (abs(a-b)>tol) then
         print *,trim(label),a,b
         error stop 'assertion failed'
      end if
   end subroutine assert_close
end program test_distributions
