program test_prediction
  use stochvoltmb
  implicit none
  type(sv_rng_state) :: rng
  type(sv_fit_result) :: fit
  type(sv_prediction) :: pred
  type(sv_prediction_summary) :: s
  real(dp), allocatable :: draws(:,:)
  real(dp) :: probs(2)
  integer :: stat

  fit%model=sv_skew_gaussian_leverage
  fit%params%sigma_y=0.2_dp; fit%params%sigma_h=0.25_dp; fit%params%phi=0.9_dp
  fit%params%alpha=-1.5_dp; fit%params%rho=-0.5_dp
  allocate(fit%theta(5),fit%theta_cov(5,5),fit%h(40),fit%h_se(40))
  call parameters_to_theta(fit%params,fit%model,fit%theta,stat)
  fit%theta_cov=0.0_dp
  fit%theta_cov(1,1)=0.0025_dp; fit%theta_cov(2,2)=0.0025_dp
  fit%theta_cov(3,3)=0.01_dp; fit%theta_cov(4,4)=0.01_dp; fit%theta_cov(5,5)=0.01_dp
  fit%h=0.0_dp; fit%h_se=0.1_dp
  call rng%seed(9911_8)
  call simulate_parameters(fit,300,rng,draws,stat)
  call check(stat==sv_ok .and. all(shape(draws)==[5,300]),'parameter simulation')
  call predict_sv(fit,8,500,rng,pred,include_parameters=.true.)
  call check(pred%status==sv_ok,'prediction status')
  call check(all(shape(pred%y)==[8,500]) .and. all(shape(pred%h)==[8,500]),'prediction dimensions')
  call check(all(finite_real(pred%y)) .and. all(pred%h_exp>0.0_dp),'prediction values')
  probs=[0.025_dp,0.975_dp]
  call summarize_prediction(pred,probs,s,include_mean=.true.)
  call check(s%status==sv_ok,'prediction summary status')
  call check(all(shape(s%y_quantiles)==[8,2]) .and. size(s%y_mean)==8,'summary dimensions')
  call check(all(s%y_quantiles(:,1)<=s%y_quantiles(:,2)),'ordered quantiles')
  print '(a)','test_prediction: PASS'
contains
  subroutine check(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//message
      error stop 1
    end if
  end subroutine check
end program test_prediction
