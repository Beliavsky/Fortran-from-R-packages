! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fMultivar authors and translation contributors.
! This file is part of fmultivar-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
program test_compatibility
  use fmultivar, only : dp, i8, dmvnorm, pmvnorm, qmvnorm, rmvnorm, &
    dmvt, pmvt, qmvt, rmvt, dmsn, pmsn, rmsn, dmst, pmst, rmst, &
    dmsc, pmsc, rmsc, dmvsnorm, pmvsnorm, rmvsnorm, dmvst, pmvst, rmvst, &
    msn_fit, mst_fit, msc_fit, mvfit, skew_fit_result, delliptical2d, &
    integrate2d, adapt, integration_result, griddata, grid_data, squarebinning, &
    hexbinning, binning_result
  implicit none
  real(dp) :: mu(1),omega(1,1),alpha(1),point(1),lower(1),upper(1)
  real(dp) :: p,se,q
  real(dp),allocatable :: x(:,:),y(:),z(:,:)
  type(skew_fit_result) :: fit
  type(integration_result) :: int_result
  type(grid_data) :: gd
  type(binning_result) :: bins_result
  logical :: ok
  integer :: failures
  failures=0;mu=[0.0_dp];omega=reshape([1.0_dp],[1,1]);alpha=[1.5_dp]
  point=[0.2_dp];lower=-huge(1.0_dp);upper=[0.3_dp]
  call check_true(dmvnorm(point,mu,omega,ok)>0.0_dp.and.ok,'dmvnorm',failures)
  call pmvnorm(lower,upper,mu,omega,p,se,30000,1001_i8,ok)
  call check_true(ok.and.p>0.0_dp.and.p<1.0_dp.and.se<0.01_dp,'pmvnorm',failures)
  q=qmvnorm(0.5_dp,mu,omega,30000)
  call check_true(abs(q)<0.06_dp,'qmvnorm',failures)
  call rmvnorm(100,mu,omega,x,1002_i8,ok);call check_true(ok.and.size(x,1)==100,'rmvnorm',failures)

  call check_true(dmvt(point,mu,omega,6.0_dp,ok)>0.0_dp.and.ok,'dmvt',failures)
  call pmvt(lower,upper,mu,omega,6.0_dp,p,se,30000,1003_i8,ok)
  call check_true(ok.and.p>0.0_dp.and.p<1.0_dp,'pmvt',failures)
  q=qmvt(0.5_dp,mu,omega,6.0_dp,30000)
  call check_true(abs(q)<0.08_dp,'qmvt',failures)
  call rmvt(100,mu,omega,6.0_dp,x,1004_i8,ok);call check_true(ok.and.size(x,1)==100,'rmvt',failures)

  call check_true(dmsn(point,mu,omega,alpha,ok)>0.0_dp.and.ok,'dmsn',failures)
  call pmsn(lower,upper,mu,omega,alpha,p,se,25000,1005_i8,ok)
  call check_true(ok.and.p>=0.0_dp.and.p<=1.0_dp,'pmsn',failures)
  call rmsn(180,mu,omega,alpha,x,1006_i8,ok);call check_true(ok,'rmsn',failures)
  fit=msn_fit(x,600,1.0e-5_dp);call check_true(fit%loglik>-huge(1.0_dp)/2.0_dp,'msn_fit',failures)

  call check_true(dmst(point,mu,omega,alpha,5.0_dp,ok)>0.0_dp.and.ok,'dmst',failures)
  call pmst(lower,upper,mu,omega,alpha,5.0_dp,p,se,25000,1007_i8,ok)
  call check_true(ok.and.p>=0.0_dp.and.p<=1.0_dp,'pmst',failures)
  call rmst(180,mu,omega,alpha,5.0_dp,x,1008_i8,ok);call check_true(ok,'rmst',failures)
  fit=mst_fit(x,5.0_dp,700,1.0e-5_dp)
  call check_true(abs(fit%nu-5.0_dp)<1.0e-12_dp,'mst_fit',failures)

  call check_true(dmsc(point,mu,omega,alpha,ok)>0.0_dp.and.ok,'dmsc',failures)
  call pmsc(lower,upper,mu,omega,alpha,p,se,25000,1009_i8,ok)
  call check_true(ok.and.p>=0.0_dp.and.p<=1.0_dp,'pmsc',failures)
  call rmsc(180,mu,omega,alpha,x,1010_i8,ok);call check_true(ok,'rmsc',failures)
  fit=msc_fit(x,700,1.0e-5_dp);call check_true(abs(fit%nu-1.0_dp)<1.0e-12_dp,'msc_fit',failures)
  fit=mvfit(x,'normal');call check_true(fit%converged,'mvfit',failures)

  call check_true(abs(dmvsnorm(point,mu,omega,alpha)-dmsn(point,mu,omega,alpha))<1.0e-14_dp, &
    'dmvsnorm alias',failures)
  call pmvsnorm(lower,upper,mu,omega,alpha,p,se,15000,1011_i8,ok)
  call check_true(ok.and.p>=0.0_dp.and.p<=1.0_dp,'pmvsnorm alias',failures)
  call rmvsnorm(40,mu,omega,alpha,x,1012_i8,ok);call check_true(ok,'rmvsnorm alias',failures)
  call check_true(abs(dmvst(point,mu,omega,alpha,5.0_dp)-dmst(point,mu,omega,alpha,5.0_dp))<1.0e-14_dp, &
    'dmvst alias',failures)
  call pmvst(lower,upper,mu,omega,alpha,5.0_dp,p,se,15000,1013_i8,ok)
  call check_true(ok.and.p>=0.0_dp.and.p<=1.0_dp,'pmvst alias',failures)
  call rmvst(40,mu,omega,alpha,5.0_dp,x,1014_i8,ok);call check_true(ok,'rmvst alias',failures)
  call check_true(delliptical2d(0.2_dp,-0.1_dp,0.3_dp,'laplace')>0.0_dp, &
    'delliptical2d alias',failures)
  int_result=integrate2d(product_fun,1.0e-7_dp)
  call check_true(abs(int_result%value-0.25_dp)<2.0e-6_dp,'integrate2d alias',failures)
  int_result=adapt([0.0_dp,0.0_dp],[1.0_dp,1.0_dp],sum_fun,2.0e-3_dp,32768)
  call check_true(abs(int_result%value-1.0_dp)<2.0e-3_dp,'adapt alias',failures)
  allocate(z(2,2));z=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[2,2])
  gd=griddata([0.0_dp,1.0_dp],[0.0_dp,1.0_dp],z)
  call check_true(all(abs(gd%z-z)<1.0e-14_dp),'griddata alias',failures)
  allocate(y(size(x,1)));y=0.5_dp*x(:,1)
  bins_result=squarebinning(x(:,1),y,8,7)
  call check_true(sum(bins_result%count)==size(x,1),'squarebinning alias',failures)
  bins_result=hexbinning(x(:,1),y,8)
  call check_true(sum(bins_result%count)==size(x,1),'hexbinning alias',failures)

  if(failures>0)then
    write(*,'(a,i0)')'Compatibility test failures: ',failures;error stop 1
  end if
  write(*,'(a)')'Original-name distribution and fitting compatibility tests passed.'
contains
  function product_fun(a,b) result(v)
    real(dp),intent(in)::a,b
    real(dp)::v
    v=a*b
  end function product_fun
  function sum_fun(a) result(v)
    real(dp),intent(in)::a(:)
    real(dp)::v
    v=sum(a)
  end function sum_fun
  subroutine check_true(cond,name,nfail)
    logical,intent(in)::cond
    character(len=*),intent(in)::name
    integer,intent(inout)::nfail
    if(.not.cond)then;write(*,'(a)')trim(name)//' failed';nfail=nfail+1;end if
  end subroutine check_true
end program test_compatibility
