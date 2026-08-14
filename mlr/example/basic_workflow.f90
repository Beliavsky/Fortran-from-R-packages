program basic_workflow
  use mlr_kinds, only : dp, i8
  use mlr_rng, only : rng_state, rng_seed
  use mlr_types, only : resample_plan, metric_summary, linear_model
  use mlr_resampling, only : make_kfold
  use mlr_evaluate, only : resample_regression
  use mlr_metrics, only : measure_rmse
  use mlr_learners, only : fit_linear_regression, predict_linear_regression
  implicit none
  real(dp)::x(50,2),y(50);integer::i
  type(rng_state)::rng;type(resample_plan)::cv;type(metric_summary)::perf
  do i=1,50
    x(i,1)=real(i,dp)/50.0_dp;x(i,2)=cos(0.2_dp*real(i,dp));y(i)=3.0_dp+1.5_dp*x(i,1)-0.75_dp*x(i,2)
  end do
  call rng_seed(rng,2026_i8);call make_kfold(50,5,rng,cv);call resample_regression(x,y,cv,predict_lm,measure_rmse,perf)
  print '(a,f10.6,a,f10.6)', '5-fold RMSE mean = ',perf%mean,', sd = ',perf%sd
contains
  subroutine predict_lm(xtrain,ytrain,xtest,pred)
    real(dp),intent(in)::xtrain(:,:),ytrain(:),xtest(:,:);real(dp),allocatable,intent(out)::pred(:)
    type(linear_model)::m
    call fit_linear_regression(xtrain,ytrain,m);call predict_linear_regression(m,xtest,pred)
  end subroutine
end program basic_workflow
