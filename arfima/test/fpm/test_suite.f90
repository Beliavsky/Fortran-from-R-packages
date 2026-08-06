program test_suite
  use arfima
  implicit none
  real(dp),allocatable::g(:),p(:)
  type(arfima_error)::err
  p=ar_to_pacf([0.3_dp,-0.1_dp])
  if(maxval(abs(pacf_to_ar(p)-[0.3_dp,-0.1_dp]))>1.0e-12_dp) error stop 1
  call tacvf_arma([0.5_dp],[real(dp)::],3,1.0_dp,g,err)
  if(err%code/=arfima_ok) error stop 2
  if(abs(g(1)-4.0_dp/3.0_dp)>1.0e-12_dp) error stop 3
  write(*,'(a)') 'arfima-test-suite: PASS'
end program test_suite
