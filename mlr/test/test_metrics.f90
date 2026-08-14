program test_metrics
  use mlr_kinds, only : dp
  use mlr_metrics
  implicit none
  real(dp) :: y(4), p(4), pr(4), prob(4,2)
  integer :: truth(4), pred(4)
  y=[1.0_dp,2.0_dp,3.0_dp,4.0_dp];p=[1.0_dp,2.0_dp,4.0_dp,3.0_dp]
  if(abs(measure_sse(y,p)-2.0_dp)>1.0e-12_dp)error stop 'sse'
  if(abs(measure_mse(y,p)-0.5_dp)>1.0e-12_dp)error stop 'mse'
  if(abs(measure_rmse(y,p)-sqrt(0.5_dp))>1.0e-12_dp)error stop 'rmse'
  if(abs(measure_mae(y,p)-0.5_dp)>1.0e-12_dp)error stop 'mae'
  if(abs(measure_rsq(y,p)-0.6_dp)>1.0e-12_dp)error stop 'rsq'
  if(measure_spearman(y,y)<0.999999999_dp)error stop 'spearman'
  if(measure_kendall(y,y)<0.999999999_dp)error stop 'kendall'
  truth=[1,1,2,2];pred=[1,2,2,2]
  if(abs(measure_acc(truth,pred)-0.75_dp)>1.0e-12_dp)error stop 'acc'
  if(abs(measure_bac(truth,pred,2)-0.75_dp)>1.0e-12_dp)error stop 'bac'
  if(abs(measure_f1(truth,pred,2)-0.8_dp)>1.0e-12_dp)error stop 'f1'
  pr=[0.1_dp,0.4_dp,0.35_dp,0.8_dp]
  truth=[1,2,1,2]
  if(abs(measure_auc(pr,truth,2)-1.0_dp)>1.0e-12_dp)error stop 'auc'
  prob(:,1)=1.0_dp-pr;prob(:,2)=pr
  if(measure_logloss(prob,truth)<=0.0_dp)error stop 'logloss'
  print *, 'test_metrics: PASS'
end program test_metrics
