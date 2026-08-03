! SPDX-License-Identifier: GPL-2.0-only
program test_mmd_misc
  use kernlab
  implicit none
  real(dp) :: x(10,1), y(10,1)
  type(kernel_spec) :: ker
  type(mmd_result) :: mr
  type(kpca_result) :: kh
  type(kernel_model) :: kf
  real(dp), allocatable :: s(:,:)
  integer :: i,st
  do i=1,10
    x(i,1)=0.1_dp*real(i,dp);y(i,1)=2.0_dp+0.1_dp*real(i,dp)
  end do
  ker=rbfdot(1.0_dp)
  call kmmd(x,y,ker,mr,bootstrap=.true.,ntimes=20)
  call check(mr%status==KL_SUCCESS.and.mr%mmd1>0.2_dp,'kmmd')
  call kha(x,ker,kh,features=2,maxiter=10)
  call check(kh%status==KL_SUCCESS,'kha')
  call kfa(x,ker,kf,features=3,subset=8)
  call kfa_predict(kf,x,s,st)
  call check(st==KL_SUCCESS.and.size(s,2)>=2,'kfa')
  print '(a)', 'test_mmd_misc: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok;character(len=*),intent(in)::msg
    if(.not.ok)then;write(*,'(a)')'FAIL: '//msg;error stop 1;end if
  end subroutine check
end program test_mmd_misc
