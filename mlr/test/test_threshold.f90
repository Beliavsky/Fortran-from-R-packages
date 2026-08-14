program test_threshold
  use mlr_kinds, only : dp
  use mlr_threshold, only : tune_binary_threshold
  use mlr_metrics, only : measure_acc
  implicit none
  real(dp)::p(6),th,perf;integer::y(6)
  p=[0.1_dp,0.2_dp,0.4_dp,0.6_dp,0.8_dp,0.9_dp];y=[1,1,1,2,2,2]
  call tune_binary_threshold(p,y,1,2,acc_cb,.false.,th,perf)
  if(abs(perf-1.0_dp)>1.0e-12_dp)error stop 'threshold performance'
  if(th<=0.4_dp.or.th>=0.6_dp)error stop 'threshold location'
  print *, 'test_threshold: PASS'
contains
  real(dp) function acc_cb(truth,response)
    integer,intent(in)::truth(:),response(:);acc_cb=measure_acc(truth,response)
  end function
end program test_threshold
