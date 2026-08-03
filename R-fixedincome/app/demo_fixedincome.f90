! SPDX-License-Identifier: MIT
program demo_fixedincome
   use fixedincome
   implicit none
   type(spot_rate_curve_t) :: curve
   type(forward_rate_t) :: forwards
   type(term_t) :: terms
   real(dp), allocatable :: rates(:), discounts(:), interpolated(:)
   integer :: i, status

   terms = term([30.0_dp, 90.0_dp, 180.0_dp, 365.0_dp, 730.0_dp], 'days')
   rates = [0.049_dp, 0.050_dp, 0.0515_dp, 0.053_dp, 0.054_dp]
   curve = spotratecurve(rates, terms, 'discrete', 'actual/365', 'actual', &
                         gregorian_to_ordinal(2026, 7, 30), status)
   call set_interpolation(curve, interp_flatforward(), status)

   discounts = curve_discount(curve, status)
   interpolated = interpolate(curve, [45.0_dp, 120.0_dp, 270.0_dp, 540.0_dp], status)
   forwards = forwardrate_from_curve(curve, status)

   print '(a)', 'Term(days)  Spot rate   Discount factor'
   do i = 1, curve%size()
      print '(f10.1,2x,f9.6,2x,f14.10)', curve%term_days(i), curve%rate(i), discounts(i)
   end do
   print '(/,a)', 'Flat-forward interpolated rates:'
   print '(*(f10.6,1x))', interpolated
   print '(/,a)', 'Successive forward rates:'
   print '(*(f10.6,1x))', forwards%rate
end program demo_fixedincome
