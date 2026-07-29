! SPDX-License-Identifier: GPL-3.0-or-later
program test_fit
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use acdm
  implicit none
  integer,parameter::n=450
  real(dp)::x(n),true_par(3),forecast(5),bootse(3)
  real(dp),allocatable::score(:,:)
  logical::fixed(3)
  type(acd_fit_result)::fixed_fit
  type(acd_order)::ord
  type(rng_state)::rng
  type(acd_fit_options)::opt
  type(acd_fit_result)::fit
  integer::st,failures
  failures=0;ord=acd_order(1,0,1);true_par=[0.20_dp,0.15_dp,0.70_dp]
  call seed_rng(rng,987654)
  call simulate_acd(n,MODEL_ACD,ord,true_par,DIST_EXPONENTIAL,[real(dp)::],x,st,rng,burn=300)
  call assert_true(st==ACDM_SUCCESS,'simulation')
  opt%model=MODEL_ACD;opt%dist=DIST_EXPONENTIAL;opt%order=ord
  opt%max_iterations=2500;opt%restarts=2;opt%tolerance=1e-8_dp;opt%seed=8123
  opt%compute_hessian=.true.;opt%compute_robust_se=.true.
  call acd_fit_model(x,opt,fit)
  call assert_true(fit%status==ACDM_SUCCESS,'fit status')
  if(fit%status==ACDM_SUCCESS)then
    call assert_true(fit%convergence==0,'fit convergence')
    call assert_true(all(ieee_is_finite(fit%parameters)),'finite parameters')
    call assert_true(fit%parameters(1)>0._dp,'positive omega')
    call assert_true(fit%parameters(2)>=-0.05_dp.and.fit%parameters(2)<0.5_dp,'alpha range')
    call assert_true(fit%parameters(3)>0.3_dp.and.fit%parameters(3)<0.95_dp,'beta range')
    call assert_true(fit%parameters(2)+fit%parameters(3)<1.15_dp,'persistence')
    call assert_true(abs(fit%parameters(1)-true_par(1))<0.18_dp,'omega recovery')
    call assert_true(abs(fit%parameters(2)-true_par(2))<0.15_dp,'alpha recovery')
    call assert_true(abs(fit%parameters(3)-true_par(3))<0.20_dp,'beta recovery')
    call assert_true(allocated(fit%hessian).and.allocated(fit%covariance),'hessian output')
    call assert_true(all(fit%standard_errors>=0._dp),'standard errors')
    call forecast_acd(fit,opt,5,forecast,st)
    call assert_true(st==ACDM_SUCCESS.and.all(forecast>0._dp),'forecast')
    call assert_true(fit%aic>0._dp.and.fit%bic>fit%aic,'information criteria')
    call acd_score_matrix(x,opt,fit%parameters,score,st)
    call assert_true(st==ACDM_SUCCESS.and.size(score,1)==n.and.size(score,2)==3,'score matrix')
    fixed=[.false.,.false.,.true.]
    call acd_fit_model(x,opt,fixed_fit,start=true_par,fixed_mask=fixed)
    call assert_true(fixed_fit%status==ACDM_SUCCESS,'fixed parameter fit')
    if(fixed_fit%status==ACDM_SUCCESS) call assert_true(abs(fixed_fit%parameters(3)-true_par(3))<1e-14_dp,'fixed beta')
    call acd_bootstrap_se(x,opt,fit,4,bootse,st)
    call assert_true(st==ACDM_SUCCESS.and.all(bootse>=0._dp),'bootstrap standard errors')
  end if
  if(failures>0)error stop 'test_fit failed'
  print '(a,3f12.6)','test_fit: PASS parameters=',fit%parameters
contains
  subroutine assert_true(ok,label)
    logical,intent(in)::ok;character(*),intent(in)::label
    if(.not.ok)then;failures=failures+1;print *,'FAIL ',trim(label);end if
  end subroutine
end program
