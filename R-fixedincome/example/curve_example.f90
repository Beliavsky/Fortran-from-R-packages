! SPDX-License-Identifier: MIT
program curve_example
   use fixedincome
   implicit none
   type(spot_rate_curve_t) :: curve
   type(forward_rate_t) :: forwards
   real(dp), allocatable :: rates(:)
   integer :: status

   curve = spotratecurve([0.0719_dp, 0.056_dp, 0.0674_dp, 0.0687_dp, 0.07_dp], &
      term([1.0_dp, 11.0_dp, 26.0_dp, 27.0_dp, 28.0_dp], 'days'), &
      'discrete', 'actual/365', 'actual', status=status)
   call set_interpolation(curve, interp_linear(), status)
   rates = interpolate(curve, [5.0_dp, 10.0_dp, 20.0_dp], status)
   forwards = forwardrate_from_curve(curve, status)

   print '(a,*(f10.6,1x))', 'interpolated: ', rates
   print '(a,*(f10.6,1x))', 'forwards:     ', forwards%rate
end program curve_example
