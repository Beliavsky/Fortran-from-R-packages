! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
program test_models_jumps
  use highfrequency, only: dp, fit_har, har_forecast, har_model
  use highfrequency, only: fit_heavy, heavy_model, heavy_forecast
  use highfrequency, only: bns_jump_test, iv_inference
  use highfrequency, only: jump_test_result, iv_inference_result
  implicit none
  real(dp) :: series(80), returns(120), rm(120), forecast, hf(3,2)
  integer :: periods(3), i
  type(har_model) :: har
  type(heavy_model) :: heavy
  type(jump_test_result) :: jump
  type(iv_inference_result) :: infer

  series(1)=1.0_dp
  do i=2,size(series)
    series(i)=0.2_dp+0.7_dp*series(i-1)+0.01_dp*sin(real(i,dp))
  end do
  periods=[1,5,10]
  har=fit_har(series,periods)
  if(.not.har%fitted_ok)error stop 1
  forecast=har_forecast(har,series,periods)
  if(forecast<=0.0_dp)error stop 1

  do i=1,size(returns)
    returns(i)=0.01_dp*sin(0.31_dp*real(i,dp))+0.004_dp*cos(0.13_dp*real(i,dp))
    rm(i)=0.00015_dp+0.00004_dp*(1.0_dp+sin(0.11_dp*real(i,dp)))
  end do
  heavy=fit_heavy(returns,rm)
  if(.not.heavy%fitted_ok)error stop 1
  if(heavy%alpha+heavy%beta>=1.0_dp)error stop 1
  hf=heavy_forecast(heavy,rm(size(rm)),3)
  if(any(hf<=0.0_dp))error stop 1

  jump=bns_jump_test(returns)
  if(jump%p_value<0.0_dp .or. jump%p_value>1.0_dp)error stop 1
  infer=iv_inference(returns)
  if(infer%standard_error<0.0_dp .or. infer%upper<infer%lower)error stop 1
  print '(a)', 'test_models_jumps: PASS'
end program test_models_jumps
