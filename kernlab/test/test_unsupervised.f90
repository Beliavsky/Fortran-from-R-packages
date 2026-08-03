! SPDX-License-Identifier: GPL-2.0-only
program test_unsupervised
  use kernlab
  implicit none
  real(dp) :: x(8,2), y(8,2), labels_real(8)
  type(kernel_spec) :: ker
  type(kpca_result) :: pc
  type(kcca_result) :: cc
  type(cluster_result) :: cl, sp
  type(ranking_result) :: rr
  type(csi_result) :: cr
  real(dp), allocatable :: pred(:,:)
  integer :: st

  x = reshape([ -2.0_dp,-2.0_dp, -1.8_dp,-2.2_dp, -2.2_dp,-1.8_dp, -1.9_dp,-1.7_dp, &
                 2.0_dp, 2.0_dp,  1.8_dp, 2.2_dp,  2.2_dp, 1.8_dp,  1.9_dp, 1.7_dp ],[8,2],order=[2,1])
  y(:,1)=x(:,1)+0.1_dp*x(:,2);y(:,2)=x(:,2)-0.1_dp*x(:,1)
  ker=rbfdot(0.5_dp)
  call kpca(x,ker,pc,features=3)
  call check(pc%status==KL_SUCCESS.and.size(pc%rotated,2)>=2,'kpca')
  call kpca_predict(pc,x,pred,st)
  call check(st==KL_SUCCESS.and.maxval(abs(pred-pc%rotated))<1.0e-7_dp,'kpca predict')
  call kcca(x,y,ker,cc,gamma=0.1_dp,ncomps=2)
  call check(cc%status==KL_SUCCESS.and.size(cc%correlations)==2,'kcca')
  call kkmeans(x,2,ker,cl,maxiter=100)
  call check(cl%status==KL_SUCCESS.and.count(cl%labels(1:4)==cl%labels(1))>=3,'kkmeans')
  call specc(x,2,ker,sp,iterations=100)
  call check(sp%status==KL_SUCCESS.and.count(sp%labels(1:4)==sp%labels(1))>=3,'specc')
  labels_real=0.0_dp;labels_real(1)=1.0_dp
  call ranking(x,labels_real,ker,rr,iterations=20)
  call check(rr%status==KL_SUCCESS.and.size(rr%score)==8,'ranking')
  call csi(x,y,ker,3,cr)
  call check(cr%status==KL_SUCCESS.and.cr%rank>=2,'csi')
  print '(a)', 'test_unsupervised: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok;character(len=*),intent(in)::msg
    if(.not.ok)then;write(*,'(a)')'FAIL: '//msg;error stop 1;end if
  end subroutine check
end program test_unsupervised
