! SPDX-License-Identifier: GPL-2.0-or-later
program test_parametric
  use mixtools
  implicit none
  type(em_control) :: ctl
  type(mixture_result) :: nfit, rfit
  type(mv_mixture_result) :: mvfit
  type(gamma_mixture_result) :: gfit
  type(multinomial_mixture_result) :: mfit
  type(rng_state) :: rng
  real(dp), allocatable :: x(:), gx(:), mvx(:,:), reps(:,:)
  real(dp) :: ycat(8,3)
  integer :: i, status

  ctl%max_iterations=500;ctl%tolerance=1.0e-7_dp;call rng_seed(rng,13579)
  allocate(x(200));call rnormmix(rng,200,[0.45_dp,0.55_dp],[-2.0_dp,2.5_dp],[0.5_dp,0.7_dp],x)
  call normalmixEM(x,2,nfit,ctl)
  call check(nfit%status==0.and.nfit%converged,"normal mixture convergence")
  call check(abs(sum(nfit%lambda)-1.0_dp)<1.0e-10_dp,"normal weights")
  call check(minval(nfit%mu)<-1.0_dp.and.maxval(nfit%mu)>1.5_dp,"normal means")

  allocate(gx(250));do i=1,125;gx(i)=random_gamma(rng,2.0_dp,0.7_dp);end do
  do i=126,250;gx(i)=random_gamma(rng,7.0_dp,0.5_dp);end do
  call gammamixEM(gx,2,gfit,ctl)
  call check(gfit%status==0.or.gfit%status==MIXTOOLS_NOT_CONVERGED,"gamma mixture status")
  call check(all(gfit%shape>0.0_dp).and.all(gfit%scale>0.0_dp),"gamma parameters")

  allocate(mvx(160,2));call rmvnormmix(rng,160,[0.5_dp,0.5_dp], &
    reshape([-2.0_dp,-2.0_dp,2.0_dp,2.0_dp],[2,2]),make_covariances(),mvx,status=status)
  call check(status==0,"multivariate simulation")
  call mvnormalmixEM(mvx,2,mvfit,ctl)
  call check(mvfit%status==0.or.mvfit%status==MIXTOOLS_NOT_CONVERGED,"MV normal mixture")
  call check(size(mvfit%sigma,3)==2,"MV covariance count")

  ycat=reshape([8.0_dp,1.0_dp,1.0_dp,7.0_dp,2.0_dp,1.0_dp,9.0_dp,1.0_dp,0.0_dp, &
    6.0_dp,3.0_dp,1.0_dp,1.0_dp,2.0_dp,7.0_dp,0.0_dp,1.0_dp,9.0_dp,2.0_dp,1.0_dp,7.0_dp,1.0_dp,1.0_dp,8.0_dp],[8,3])
  call multmixEM(ycat,2,mfit,ctl)
  call check(mfit%status==0.and.maxval(abs(sum(mfit%theta,dim=2)-1.0_dp))<1.0e-10_dp,"multinomial mixture")

  allocate(reps(3,120));do i=1,60;reps(:,i)=-1.5_dp+[random_normal(rng),random_normal(rng),random_normal(rng)]*0.4_dp;end do
  do i=61,120;reps(:,i)=2.0_dp+[random_normal(rng),random_normal(rng),random_normal(rng)]*0.5_dp;end do
  call repnormmixEM(reps,2,rfit,ctl)
  call check(rfit%status==0.and.minval(rfit%mu)<0.0_dp.and.maxval(rfit%mu)>1.0_dp,"repeated normal mixture")
  print '(a)', 'test_parametric: PASS'
contains
  function make_covariances() result(s)
    real(dp)::s(2,2,2)
    s=0.0_dp;s(1,1,:)=0.5_dp;s(2,2,:)=0.5_dp;s(1,2,:)=0.1_dp;s(2,1,:)=0.1_dp
  end function make_covariances
  subroutine check(condition,message)
    logical,intent(in)::condition;character(len=*),intent(in)::message
    if(.not.condition)then;write(*,'(a)')'FAIL: '//message;error stop 1;end if
  end subroutine check
end program test_parametric
