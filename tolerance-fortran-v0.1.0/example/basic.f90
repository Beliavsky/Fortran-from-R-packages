program basic
  use tolerance
  implicit none

  real(dp) :: x(10)
  type(tolerance_interval) :: normal_ti, np_ti

  x = [-1.2_dp, -0.8_dp, -0.4_dp, -0.1_dp, 0.0_dp, &
        0.2_dp,  0.5_dp,  0.8_dp,  1.1_dp, 1.5_dp]

  normal_ti = normtol_int(x, alpha=0.05_dp, p=0.90_dp, side=2, method='HE')
  np_ti = nptol_int(x, alpha=0.10_dp, p=0.80_dp, side=1, method='WILKS')

  print '(a,2f12.6)', 'Normal 2-sided limits: ', normal_ti%lower, normal_ti%upper
  print '(a,2f12.6)', 'Wilks 1-sided limits: ', np_ti%lower, np_ti%upper
end program basic
