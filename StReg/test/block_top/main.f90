! SPDX-License-Identifier: GPL-2.0-only
program test_block_top
   use streg, only : dp, block_top_covariance, streg_success
   use test_support_block_top
   implicit none
   real(dp) :: a(4),expected11,expected12
   real(dp), allocatable :: s(:,:)
   integer :: status
   a=[1.0_dp,2.0_dp,3.0_dp,4.0_dp]
   call block_top_covariance(a,1,1,s,status)
   expected11=1.0_dp+4.0_dp+9.0_dp+16.0_dp+1.0_dp
   expected12=2.0_dp+6.0_dp+12.0_dp+4.0_dp+1.0_dp
   call assert_true(status==streg_success,'block covariance status')
   call assert_true(all(shape(s)==[2,2]),'block covariance dimensions')
   call assert_close(s(1,1),expected11,1.0e-13_dp,'first diagonal fixture')
   call assert_close(s(2,2),expected11,1.0e-13_dp,'second diagonal fixture')
   call assert_close(s(1,2),expected12,1.0e-13_dp,'off-diagonal fixture')
   call assert_close(s(2,1),expected12,1.0e-13_dp,'symmetry fixture')
   write(*,'(a)')'test_block_top: PASS'
end program test_block_top
