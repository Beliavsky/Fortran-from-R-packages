program example_smith_wilson
   use smith_wilson_yield_curve, only : dp, fit_smith_wilson_curve_to_instruments, &
                                        make_market_instrument, market_instrument, &
                                        smith_wilson_curve, sw_success
   implicit none

   type(market_instrument) :: instruments(4)
   type(smith_wilson_curve) :: curve
   real(dp) :: terms(8)
   real(dp), allocatable :: discounts(:), spot_rates(:), repriced(:)
   integer :: i, info
   character(len=256) :: message

   call make_market_instrument('LIBOR', 1.0_dp, 0.010_dp, instruments(1), info, message)
   call make_market_instrument('SWAP', 2.0_dp, 0.020_dp, instruments(2), info, message, &
                               frequency=1.0_dp)
   call make_market_instrument('SWAP', 3.0_dp, 0.026_dp, instruments(3), info, message, &
                               frequency=1.0_dp)
   call make_market_instrument('SWAP', 5.0_dp, 0.034_dp, instruments(4), info, message, &
                               frequency=1.0_dp)

   call fit_smith_wilson_curve_to_instruments(instruments, 0.042_dp, 0.1_dp, curve, &
                                               info, message)
   if (info /= sw_success) then
      print '(a)', 'fit failed: '//trim(message)
      error stop 1
   end if

   terms = [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, 10.0_dp, 20.0_dp, 50.0_dp, 100.0_dp]
   discounts = curve%discount(terms)
   spot_rates = curve%continuous_spot(terms)
   repriced = curve%repriced_values()

   print '(a)', ' Smith-Wilson curve'
   print '(a)', '       term          discount       continuous spot'
   do i = 1, size(terms)
      print '(3f18.10)', terms(i), discounts(i), spot_rates(i)
   end do

   print '(/,a,*(f14.10,1x))', 'Repriced instruments: ', repriced
end program example_smith_wilson
