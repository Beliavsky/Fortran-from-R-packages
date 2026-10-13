program credit_curve
  use creditrisk, only: dp, credit_curve_t, at1p, cds, credit_ok
  use creditrisk_calibration, only: calibration_result_t, calibrate_cds
  implicit none
  real(dp), parameter :: time(4) = [0.5_dp, 1.0_dp, 2.0_dp, 4.0_dp]
  real(dp), parameter :: rate(4) = [0.015_dp, 0.017_dp, 0.022_dp, 0.028_dp]
  type(credit_curve_t) :: structural, swap
  type(calibration_result_t) :: fit
  integer :: i
  structural = at1p(1.0_dp, 0.65_dp, 0.4_dp, [0.19_dp, 0.23_dp, 0.21_dp, 0.26_dp], rate, time)
  if (structural%status /= credit_ok) error stop 'AT1P failed'
  swap = cds(time, structural%values(:,5), rate)
  if (swap%status /= credit_ok) error stop 'CDS failed'
  print '(a)', 'Years       Survival       CDS decimal rate'
  do i = 1, size(time)
    print '(f6.2,2f18.10)', time(i), structural%values(i,4), swap%values(i,5)
  end do
  fit = calibrate_cds(rate, time, [0.75_dp, 1.5_dp, 3.0_dp], [0.008_dp, 0.013_dp, 0.02_dp])
  if (fit%status /= credit_ok .or. .not. fit%optimizer%success) error stop 'CDS calibration failed'
  print '(a,3f14.8)', 'Calibrated hazard knots: ', fit%parameters
  print '(a,es14.5)', 'Mean squared spread error: ', fit%error
end program credit_curve
