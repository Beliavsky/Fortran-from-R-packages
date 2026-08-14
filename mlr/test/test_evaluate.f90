program test_evaluate
  use mlr_kinds, only : dp, i8
  use mlr_rng, only : rng_state, rng_seed
  use mlr_types, only : resample_plan, metric_summary, linear_model
  use mlr_resampling, only : make_kfold
  use mlr_evaluate, only : resample_regression
  use mlr_metrics, only : measure_rmse
  use mlr_learners, only : fit_linear_regression, predict_linear_regression
  implicit none
  real(dp)::x(30,2),y(30);integer::i
  type(rng_state)::rng;type(resample_plan)::plan;type(metric_summary)::res
  do i=1,30;x(i,1)=real(i,dp)/30.0_dp;x(i,2)=sin(real(i,dp));y(i)=1.5_dp+2.0_dp*x(i,1)-0.3_dp*x(i,2);end do
  call rng_seed(rng,42_i8);call make_kfold(30,5,rng,plan);call resample_regression(x,y,plan,predictor,measure_rmse,res)
  if(res%mean>1.0e-8_dp)error stop 'cv regression'
  print *, 'test_evaluate: PASS'
contains
  subroutine predictor(xtrain,ytrain,xtest,pred)
    real(dp),intent(in)::xtrain(:,:),ytrain(:),xtest(:,:);real(dp),allocatable,intent(out)::pred(:)
    type(linear_model)::m
    call fit_linear_regression(xtrain,ytrain,m);call predict_linear_regression(m,xtest,pred)
  end subroutine
end program test_evaluate
