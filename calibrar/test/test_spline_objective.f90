program test_spline_objective
  use calibrar, only : dp, spline_result, spline_par, calibration_term, calibration_objective_value
  implicit none
  type(spline_result) :: sp
  type(calibration_term) :: terms(2)
  real(dp) :: obs(4),sim(4),vals(2),tot
  call spline_par([0.0_dp,1.0_dp,2.0_dp],4,sp)
  if(maxval(abs(sp%x-[0.25_dp,0.75_dp,1.25_dp,1.75_dp]))>1.0e-10_dp) error stop "spline failed"
  obs=[1.0_dp,2.0_dp,3.0_dp,4.0_dp];sim=[1.0_dp,1.0_dp,2.0_dp,4.0_dp]
  terms(1)%first=1;terms(1)%last=2;terms(1)%fitness_type="norm2";terms(1)%weight=2.0_dp
  terms(2)%first=3;terms(2)%last=4;terms(2)%fitness_type="norm2";terms(2)%weight=1.0_dp
  call calibration_objective_value(obs,sim,terms,vals,total=tot)
  if(abs(vals(1)-1.0_dp)>1.0e-12_dp .or. abs(vals(2)-1.0_dp)>1.0e-12_dp) error stop "objective terms failed"
  if(abs(tot-3.0_dp)>1.0e-12_dp) error stop "objective aggregate failed"
  print *, "PASS test_spline_objective"
end program test_spline_objective
