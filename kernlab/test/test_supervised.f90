! SPDX-License-Identifier: GPL-2.0-only
program test_supervised
  use kernlab
  implicit none
  real(dp) :: x(8,2), yr(8), xnew(2)
  integer :: yc(8), st
  type(kernel_spec) :: ker
  type(kernel_model) :: lsreg, lscls, svm, gp, rv, qr, online
  real(dp), allocatable :: p(:,:), v(:)
  integer, allocatable :: cp(:)

  x = reshape([ -2.0_dp,-2.0_dp, -1.8_dp,-2.2_dp, -2.2_dp,-1.8_dp, -1.9_dp,-1.7_dp, &
                 2.0_dp, 2.0_dp,  1.8_dp, 2.2_dp,  2.2_dp, 1.8_dp,  1.9_dp, 1.7_dp ],[8,2],order=[2,1])
  yc=[-1,-1,-1,-1,1,1,1,1];yr=2.0_dp*x(:,1)-x(:,2)
  ker=rbfdot(0.5_dp)
  call lssvm(x,yr,ker,lsreg,tau=1.0e-4_dp);call predict_kernel_model(lsreg,x,p,st)
  call check(st==KL_SUCCESS.and.sum((p(:,1)-yr)**2)/8.0_dp<0.05_dp,'lssvm regression')
  call lssvm(x,yc,ker,lscls,tau=1.0e-3_dp);call predict_kernel_model(lscls,x,p,st,cp)
  call check(count(cp==yc)>=7,'lssvm classification')
  call ksvm(x,yc,ker,svm,cost=10.0_dp,maxiter=2000);call predict_kernel_model(svm,x,p,st,cp)
  call check(count(cp==yc)>=7,'ksvm classification')
  call gausspr(x,yr,ker,gp,var=1.0e-4_dp);call gausspr_predict_variance(gp,x,p,v,st)
  call check(st==KL_SUCCESS.and.all(v>=0.0_dp),'gaussian process')
  call rvm(x,yr,ker,rv,maxiter=100);call predict_kernel_model(rv,x,p,st)
  call check(st==KL_SUCCESS.and.all(abs(p(:,1))<100.0_dp),'rvm')
  call kqr(x,yr,ker,qr,tau=0.5_dp,maxiter=500);call predict_kernel_model(qr,x,p,st)
  call check(st==KL_SUCCESS.and.all(abs(p(:,1))<100.0_dp),'kqr')
  call onlearn(x,yr,ker,online,lambda=0.5_dp,buffer_size=6)
  xnew=[0.0_dp,0.0_dp];call inlearn(online,xnew,0.0_dp,buffer_size=6)
  call predict_kernel_model(online,reshape(xnew,[1,2]),p,st)
  call check(st==KL_SUCCESS.and.size(online%train,1)==6,'online learning')
  print '(a)', 'test_supervised: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok;character(len=*),intent(in)::msg
    if(.not.ok)then;write(*,'(a)')'FAIL: '//msg;error stop 1;end if
  end subroutine check
end program test_supervised
