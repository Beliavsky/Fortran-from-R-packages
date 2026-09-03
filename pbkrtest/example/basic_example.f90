program basic_example
   use pbkrtest, only : bootstrap_p_values, bootstrap_result_t, dp, likelihood_ratio_test, &
      lrt_result_t, pbkr_success
   implicit none
   type(bootstrap_result_t) :: boot
   type(lrt_result_t) :: lrt
   real(dp) :: reference(6)
   integer :: status

   call likelihood_ratio_test(-120.0_dp, -123.0_dp, 8, 6, lrt, status)
   if (status /= pbkr_success) error stop 'likelihood-ratio calculation failed'
   print '(a,f10.5,a,i0,a,f10.6)', 'LRT = ', lrt%statistic, ', df = ', lrt%df, ', p = ', lrt%p_value

   reference = [0.8_dp, 1.7_dp, 2.6_dp, 4.0_dp, 5.2_dp, 7.1_dp]
   call bootstrap_p_values(lrt%statistic, lrt%df, reference, boot, status)
   if (status /= pbkr_success) error stop 'bootstrap calibration failed'
   print '(a,f10.6)', 'Parametric-bootstrap p = ', boot%p_bootstrap
   print '(a,f10.6)', 'Gamma moment-match p = ', boot%p_gamma
end program basic_example
