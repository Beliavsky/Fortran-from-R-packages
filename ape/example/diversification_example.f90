program diversification_example
   use ape, only : dp, diversification_time_result, diversification_time
   implicit none

   type(diversification_time_result) :: fit
   integer :: info

   call diversification_time([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], fit, info, &
      census=[1, 1, 1, 0, 1, 0], breakpoint=3.5_dp)
   if (info /= 0) error stop 'diversification_time failed'

   print '(a,f10.6)', 'constant rate: ', fit%constant_rate
   print '(a,f10.6)', 'Weibull shape: ', fit%weibull_shape
   print '(a,f10.6)', 'early rate:    ', fit%early_rate
   print '(a,f10.6)', 'late rate:     ', fit%late_rate
end program diversification_example
