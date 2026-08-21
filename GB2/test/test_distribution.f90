program test_distribution
  use gb2, only : dp,dgb2,pgb2,qgb2,moment_gb2,expected_log_gb2,variance_log_gb2, &
    skewness_log_gb2,kurtosis_log_gb2,gini_gb2
  implicit none
  integer :: fails
  fails=0
  call check(dgb2(3.1_dp,2.3_dp,4.2_dp,1.7_dp,3.4_dp),0.29775404560015145_dp,2e-13_dp,'density',fails)
  call check(pgb2(3.1_dp,2.3_dp,4.2_dp,1.7_dp,3.4_dp),0.5403417178701424_dp,2e-13_dp,'cdf',fails)
  call check(qgb2(0.1_dp,2.3_dp,4.2_dp,1.7_dp,3.4_dp),1.597651643762604_dp,2e-12_dp,'q10',fails)
  call check(qgb2(0.5_dp,2.3_dp,4.2_dp,1.7_dp,3.4_dp),2.9671458698315933_dp,2e-12_dp,'q50',fails)
  call check(qgb2(0.9_dp,2.3_dp,4.2_dp,1.7_dp,3.4_dp),5.052913032037388_dp,2e-12_dp,'q90',fails)
  call check(moment_gb2(1._dp,2.3_dp,4.2_dp,1.7_dp,3.4_dp),3.1980398622384523_dp,2e-12_dp,'mean',fails)
  call check(moment_gb2(2._dp,2.3_dp,4.2_dp,1.7_dp,3.4_dp),12.363264742795296_dp,2e-11_dp,'second moment',fails)
  call check(expected_log_gb2(2.3_dp,4.2_dp,1.7_dp,3.4_dp),1.060728020292031_dp,2e-12_dp,'Elog',fails)
  call check(variance_log_gb2(2.3_dp,1.7_dp,3.4_dp),0.2145130865689663_dp,2e-12_dp,'Vlog',fails)
  call check(skewness_log_gb2(1.7_dp,3.4_dp),-0.40411019158604056_dp,2e-11_dp,'Slog',fails)
  call check(kurtosis_log_gb2(1.7_dp,3.4_dp),0.7510425404051441_dp,2e-11_dp,'Klog',fails)
  call check(gini_gb2(2.3_dp,1.7_dp,3.4_dp),0.2446231450316234_dp,2e-8_dp,'Gini',fails)
  if(fails>0) error stop 1
  print '(a)', 'test_distribution: PASS'
contains
  subroutine check(x,y,tol,name,fails)
    real(dp),intent(in)::x,y,tol
    character(len=*),intent(in)::name
    integer,intent(inout)::fails
    if(abs(x-y)>tol*max(1._dp,abs(y))) then
    print *,trim(name),x,y
    fails=fails+1
    end if
  end subroutine
end program
