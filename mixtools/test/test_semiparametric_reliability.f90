! SPDX-License-Identifier: GPL-2.0-or-later
program test_semiparametric_reliability
  use mixtools
  implicit none
  type(em_control) :: ctl
  type(semiparametric_result) :: npfit, symfit
  type(reliability_mixture_result) :: efit, wfit
  type(model_selection_result) :: selectfit
  type(bootstrap_result) :: boot
  type(mcmc_result) :: chain
  type(rng_state) :: rng
  real(dp), allocatable :: x(:,:), v(:), time(:), xr(:,:), yr(:)
  integer, allocatable :: event(:)
  integer :: i, n

  ctl%max_iterations=100;ctl%tolerance=1.0e-5_dp;call rng_seed(rng,97531);n=100
  allocate(x(n,2),v(n));do i=1,n/2;x(i,:)=[-2.0_dp+0.4_dp*random_normal(rng),-1.0_dp+0.5_dp*random_normal(rng)]
    v(i)=-2.0_dp+0.5_dp*random_normal(rng);end do
  do i=n/2+1,n;x(i,:)=[2.0_dp+0.5_dp*random_normal(rng),1.0_dp+0.4_dp*random_normal(rng)]
    v(i)=2.0_dp+0.6_dp*random_normal(rng);end do
  call npEM(x,2,npfit,ctl,ngrid=80)
  call check((npfit%status==0.or.npfit%status==MIXTOOLS_NOT_CONVERGED).and.size(npfit%density,1)==2,"npEM")
  call spEMsymlocN01(v,symfit,ctl,location=2.0_dp,ngrid=80)
  call check(size(symfit%density,1)==2,"symmetric semiparametric mixture")

  allocate(time(160),event(160));do i=1,80;time(i)=random_exponential(rng,1.5_dp);event(i)=1;end do
  do i=81,160;time(i)=random_weibull(rng,2.0_dp,2.0_dp);event(i)=merge(1,0,time(i)<2.5_dp);time(i)=min(time(i),2.5_dp);end do
  call expRMM_EM(time,event,2,efit,ctl)
  call check(efit%status==0.and.all(efit%scale>0.0_dp),"exponential reliability mixture")
  call weibullRMM_SEM(time,event,2,wfit,ctl)
  call check((wfit%status==0.or.wfit%status==MIXTOOLS_NOT_CONVERGED).and.all(wfit%shape>0.0_dp),"Weibull reliability mixture")

  call normalmix_model_selection(v,[1,2,3],selectfit,ctl)
  call check(selectfit%best_bic>=1.and.selectfit%best_bic<=3,"model selection")
  call boot_comp(v,1,2,4,boot,ctl,123)
  call check(boot%successful==4.and.boot%p_value>0.0_dp,"bootstrap comparison")

  allocate(xr(60,1),yr(60));do i=1,60;xr(i,1)=real(i-30,dp)/15.0_dp;yr(i)=1.0_dp+0.7_dp*xr(i,1)+0.3_dp*random_normal(rng);end do
  call regmixMH(yr,xr,1,8,chain,seed=777,burnin=8,thin=1)
  call check(chain%status==0.and.size(chain%beta_draws,1)==8,"regression mixture MCMC")
  print '(a)', 'test_semiparametric_reliability: PASS'
contains
  subroutine check(condition,message)
    logical,intent(in)::condition;character(len=*),intent(in)::message
    if(.not.condition)then;write(*,'(a)')'FAIL: '//message;error stop 1;end if
  end subroutine check
end program test_semiparametric_reliability
