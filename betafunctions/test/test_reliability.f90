program test_reliability
  use betafunctions
  implicit none
  real(dp) :: x(6,3), alpha
  type(omega_result) :: om

  x = reshape([ &
    1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp, &
    1.5_dp,2.4_dp,3.2_dp,4.3_dp,5.1_dp,6.2_dp, &
    0.8_dp,1.9_dp,2.8_dp,4.1_dp,4.9_dp,6.1_dp ], [6,3])
  alpha = cronbach_alpha(x)
  if (.not. (alpha > 0.9_dp .and. alpha <= 1.01_dp)) error stop 1
  call mcdonald_omega(x, om)
  if (.not. (om%omega > 0.9_dp .and. om%omega <= 1.01_dp)) error stop 1
  if (maxval(abs(om%observed - transpose(om%observed))) > 1.0e-14_dp) error stop 1

  print '(a)', 'test_reliability: PASS'
end program test_reliability
