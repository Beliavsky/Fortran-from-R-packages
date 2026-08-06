program test_psi_location
  use robstattm, only : dp, location_scale_result, rho_value, rho_prime, rho_second, &
    scale_m, inverse_robust_r_squared, loc_scale_m
  implicit none
  real(dp) :: x(9), s, r2
  type(location_scale_result) :: fit

  x = [-2.0_dp, -1.0_dp, -0.5_dp, -0.1_dp, 0.0_dp, 0.2_dp, 0.6_dp, 1.1_dp, 25.0_dp]
  call assert_true(abs(rho_value(0.0_dp, 'bisquare', 4.685_dp)) < 1.0e-14_dp, 'rho(0)')
  call assert_true(rho_value(100.0_dp, 'bisquare', 4.685_dp) > 0.999_dp, 'bounded rho')
  call assert_true(abs(rho_prime(0.0_dp, 'huber', 1.345_dp)) < 1.0e-12_dp, 'huber derivative')
  call assert_true(rho_second(0.0_dp, 'bisquare', 4.685_dp) > 0.0_dp, 'positive curvature')
  s = scale_m(x, 0.5_dp, 'bisquare', 1.54764_dp)
  call assert_true(s > 0.0_dp .and. s < 5.0_dp, 'robust M-scale')
  call loc_scale_m(x, fit, family='bisquare', efficiency=0.85_dp)
  call assert_true(abs(fit%location) < 1.0_dp, 'robust location')
  call assert_true(fit%scale > 0.0_dp, 'location scale')
  r2 = inverse_robust_r_squared(0.25_dp, 'bisquare', 3.4434_dp)
  call assert_true(r2 >= 0.0_dp .and. r2 <= 1.0_dp, 'inverse robust R2')
  print '(a)', 'test_psi_location: PASS'
contains
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // trim(message)
      error stop 1
    end if
  end subroutine assert_true
end program test_psi_location
