! SPDX-License-Identifier: MIT
program test_prediction
   use vasicekfit, only : dp, vasicek_fit_result, prediction_result, fit_vasicek, &
      predict_link, predict_response, predict_quantiles
   implicit none
   real(dp), parameter :: y(20) = [ &
      0.009600134840871976_dp, 0.029502890965674730_dp, 0.037198232908477340_dp, &
      0.020785442980663090_dp, 0.050685547743244040_dp, 0.017640200610933460_dp, &
      0.032230717502823025_dp, 0.072262483067218770_dp, 0.016096229432628024_dp, &
      0.048090550856388070_dp, 0.031519799988685336_dp, 0.065329288703111180_dp, &
      0.015700650584264047_dp, 0.060401782241646530_dp, 0.039279530963784316_dp, &
      0.033014580456703310_dp, 0.081703730226032780_dp, 0.045890279005243326_dp, &
      0.071152806570663060_dp, 0.033199582822105970_dp ]
   real(dp) :: x(20,2), new_x(3,2), alpha(2)
   integer :: i
   type(vasicek_fit_result) :: fit
   type(prediction_result) :: link, response, median, tails

   do i = 1, 20
      x(i,1) = -1.5_dp + 3.0_dp * real(i-1,dp) / 19.0_dp
      x(i,2) = (-1.0_dp)**(i-1) * (0.2_dp + 0.03_dp * real(i-1,dp))
   end do
   fit = fit_vasicek(y, x)
   call assert_true(fit%ok)

   new_x = reshape([-1.0_dp, 0.0_dp, 1.0_dp, 0.25_dp, -0.4_dp, 0.55_dp], [3,2])
   link = predict_link(fit, new_x)
   response = predict_response(fit, new_x)
   median = predict_quantiles(fit, [0.5_dp], new_x)
   alpha = [0.1_dp, 0.9_dp]
   tails = predict_quantiles(fit, alpha, new_x)

   call assert_true(link%ok .and. response%ok .and. median%ok .and. tails%ok)
   call assert_array(link%values(:,1), [ &
      -1.9362883063191030_dp, -1.7883325440089550_dp, -1.6811270483411893_dp], 3.0e-12_dp)
   call assert_array(response%values(:,1), [ &
      0.02941011101039863_dp, 0.04047925513379026_dp, 0.05044703686071087_dp], 3.0e-12_dp)
   call assert_array(median%values(:,1), [ &
      0.02641620009340112_dp, 0.03686118615324616_dp, 0.04636911980131028_dp], 3.0e-12_dp)
   call assert_array(tails%values(:,1), [ &
      0.01310331534629394_dp, 0.01898491667467119_dp, 0.02453547189463528_dp], 3.0e-12_dp)
   call assert_array(tails%values(:,2), [ &
      0.04952870029737790_dp, 0.06661510129652340_dp, 0.08161663162574440_dp], 3.0e-12_dp)
   call assert_true(all(tails%values(:,2) > tails%values(:,1)))

   print '(a)', 'test_prediction: PASS'

contains

   subroutine assert_array(actual, expected, tolerance)
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      if (size(actual) /= size(expected)) error stop 1
      if (maxval(abs(actual - expected)) > tolerance + tolerance * maxval(abs(expected))) then
         print '(a,es25.16)', 'maximum mismatch: ', maxval(abs(actual - expected))
         error stop 1
      end if
   end subroutine assert_array

   subroutine assert_true(condition)
      logical, intent(in) :: condition
      if (.not. condition) error stop 1
   end subroutine assert_true

end program test_prediction
