! SPDX-License-Identifier: MIT
! SPDX-FileComment: Pure curve helpers corresponding to common R stats self-start models.
module r_stats_self_start
   use r_kinds, only: dp
   implicit none
   private

   public :: ss_asymp, ss_asymp_jacobian, ss_asymp_model
   public :: ss_gompertz, ss_gompertz_jacobian, ss_gompertz_model
   public :: ss_logis, ss_logis_jacobian, ss_logis_model
   public :: ss_micmen, ss_micmen_jacobian, ss_micmen_model
   public :: ss_weibull, ss_weibull_jacobian, ss_weibull_model

contains

   pure elemental real(dp) function ss_asymp(x, asymptote, initial_response, log_rate) result(value)
      !! Evaluates R's `SSasymp` asymptotic regression curve.
      real(dp), intent(in) :: x !! Predictor value.
      real(dp), intent(in) :: asymptote !! Limiting response.
      real(dp), intent(in) :: initial_response !! Response at zero.
      real(dp), intent(in) :: log_rate !! Natural logarithm of the positive rate.

      value = asymptote + (initial_response - asymptote)*exp(-exp(log_rate)*x)
   end function ss_asymp

   pure subroutine ss_asymp_model(x, parameters, values)
      !! Evaluates `SSasymp` through the callback accepted by `nls_fit`.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix whose first column is used.
      real(dp), intent(in) :: parameters(:) !! `[Asym,R0,lrc]` parameters.
      real(dp), intent(out) :: values(size(x, 1)) !! Curve values.

      values = ss_asymp(x(:, 1), parameters(1), parameters(2), parameters(3))
   end subroutine ss_asymp_model

   pure subroutine ss_asymp_jacobian(x, parameters, jacobian)
      !! Evaluates the analytic `SSasymp` Jacobian.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix whose first column is used.
      real(dp), intent(in) :: parameters(:) !! `[Asym,R0,lrc]` parameters.
      real(dp), intent(out) :: jacobian(size(x, 1), size(parameters)) !! Model Jacobian.
      real(dp) :: decay(size(x, 1)), rate

      rate = exp(parameters(3))
      decay = exp(-rate*x(:, 1))
      jacobian(:, 1) = 1.0_dp - decay
      jacobian(:, 2) = decay
      jacobian(:, 3) = -(parameters(2) - parameters(1))*decay*rate*x(:, 1)
   end subroutine ss_asymp_jacobian

   pure elemental real(dp) function ss_gompertz(x, asymptote, b2, b3) result(value)
      !! Evaluates R's `SSgompertz` growth curve.
      real(dp), intent(in) :: x !! Predictor value.
      real(dp), intent(in) :: asymptote !! Limiting response.
      real(dp), intent(in) :: b2 !! Positive displacement parameter.
      real(dp), intent(in) :: b3 !! Positive rate base.

      value = asymptote*exp(-b2*b3**x)
   end function ss_gompertz

   pure subroutine ss_gompertz_model(x, parameters, values)
      !! Evaluates `SSgompertz` through the callback accepted by `nls_fit`.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix whose first column is used.
      real(dp), intent(in) :: parameters(:) !! `[Asym,b2,b3]` parameters.
      real(dp), intent(out) :: values(size(x, 1)) !! Curve values.

      values = ss_gompertz(x(:, 1), parameters(1), parameters(2), parameters(3))
   end subroutine ss_gompertz_model

   pure subroutine ss_gompertz_jacobian(x, parameters, jacobian)
      !! Evaluates the analytic `SSgompertz` Jacobian.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix whose first column is used.
      real(dp), intent(in) :: parameters(:) !! `[Asym,b2,b3]` parameters.
      real(dp), intent(out) :: jacobian(size(x, 1), size(parameters)) !! Model Jacobian.
      real(dp) :: curve(size(x, 1)), power(size(x, 1))

      power = parameters(3)**x(:, 1)
      curve = exp(-parameters(2)*power)
      jacobian(:, 1) = curve
      jacobian(:, 2) = -parameters(1)*curve*power
      jacobian(:, 3) = -parameters(1)*curve*parameters(2)* &
                       parameters(3)**(x(:, 1) - 1.0_dp)*x(:, 1)
   end subroutine ss_gompertz_jacobian

   pure elemental real(dp) function ss_logis(x, asymptote, midpoint, scale) result(value)
      !! Evaluates R's `SSlogis` logistic growth curve.
      real(dp), intent(in) :: x !! Predictor value.
      real(dp), intent(in) :: asymptote !! Limiting response.
      real(dp), intent(in) :: midpoint !! Predictor at half the asymptote.
      real(dp), intent(in) :: scale !! Logistic scale parameter.

      value = asymptote/(1.0_dp + exp((midpoint - x)/scale))
   end function ss_logis

   pure subroutine ss_logis_model(x, parameters, values)
      !! Evaluates `SSlogis` through the callback accepted by `nls_fit`.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix whose first column is used.
      real(dp), intent(in) :: parameters(:) !! `[Asym,xmid,scal]` parameters.
      real(dp), intent(out) :: values(size(x, 1)) !! Curve values.

      values = ss_logis(x(:, 1), parameters(1), parameters(2), parameters(3))
   end subroutine ss_logis_model

   pure subroutine ss_logis_jacobian(x, parameters, jacobian)
      !! Evaluates the analytic `SSlogis` Jacobian.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix whose first column is used.
      real(dp), intent(in) :: parameters(:) !! `[Asym,xmid,scal]` parameters.
      real(dp), intent(out) :: jacobian(size(x, 1), size(parameters)) !! Model Jacobian.
      real(dp) :: exponent(size(x, 1)), denominator(size(x, 1)), midpoint_derivative(size(x, 1))

      exponent = exp((parameters(2) - x(:, 1))/parameters(3))
      denominator = 1.0_dp + exponent
      midpoint_derivative = -parameters(1)*exponent/parameters(3)/denominator**2
      jacobian(:, 1) = 1.0_dp/denominator
      jacobian(:, 2) = midpoint_derivative
      jacobian(:, 3) = midpoint_derivative*(parameters(2) - x(:, 1))/parameters(3)
   end subroutine ss_logis_jacobian

   pure elemental real(dp) function ss_micmen(x, maximum, half_saturation) result(value)
      !! Evaluates R's `SSmicmen` Michaelis-Menten curve.
      real(dp), intent(in) :: x !! Predictor value.
      real(dp), intent(in) :: maximum !! Limiting response `Vm`.
      real(dp), intent(in) :: half_saturation !! Half-saturation parameter `K`.

      value = maximum*x/(half_saturation + x)
   end function ss_micmen

   pure subroutine ss_micmen_model(x, parameters, values)
      !! Evaluates `SSmicmen` through the callback accepted by `nls_fit`.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix whose first column is used.
      real(dp), intent(in) :: parameters(:) !! `[Vm,K]` parameters.
      real(dp), intent(out) :: values(size(x, 1)) !! Curve values.

      values = ss_micmen(x(:, 1), parameters(1), parameters(2))
   end subroutine ss_micmen_model

   pure subroutine ss_micmen_jacobian(x, parameters, jacobian)
      !! Evaluates the analytic `SSmicmen` Jacobian.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix whose first column is used.
      real(dp), intent(in) :: parameters(:) !! `[Vm,K]` parameters.
      real(dp), intent(out) :: jacobian(size(x, 1), size(parameters)) !! Model Jacobian.
      real(dp) :: denominator(size(x, 1))

      denominator = parameters(2) + x(:, 1)
      jacobian(:, 1) = x(:, 1)/denominator
      jacobian(:, 2) = -parameters(1)*x(:, 1)/denominator**2
   end subroutine ss_micmen_jacobian

   pure elemental real(dp) function ss_weibull(x, asymptote, drop, log_rate, power) result(value)
      !! Evaluates R's `SSweibull` growth curve.
      real(dp), intent(in) :: x !! Nonnegative predictor value.
      real(dp), intent(in) :: asymptote !! Limiting response.
      real(dp), intent(in) :: drop !! Difference from the asymptote at zero.
      real(dp), intent(in) :: log_rate !! Natural logarithm of the positive rate.
      real(dp), intent(in) :: power !! Positive predictor exponent.

      value = asymptote - drop*exp(-exp(log_rate)*x**power)
   end function ss_weibull

   pure subroutine ss_weibull_model(x, parameters, values)
      !! Evaluates `SSweibull` through the callback accepted by `nls_fit`.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix whose first column is used.
      real(dp), intent(in) :: parameters(:) !! `[Asym,Drop,lrc,pwr]` parameters.
      real(dp), intent(out) :: values(size(x, 1)) !! Curve values.

      values = ss_weibull(x(:, 1), parameters(1), parameters(2), parameters(3), parameters(4))
   end subroutine ss_weibull_model

   pure subroutine ss_weibull_jacobian(x, parameters, jacobian)
      !! Evaluates the analytic `SSweibull` Jacobian for positive predictors.
      real(dp), intent(in) :: x(:, :) !! Positive predictor values in the first column.
      real(dp), intent(in) :: parameters(:) !! `[Asym,Drop,lrc,pwr]` parameters.
      real(dp), intent(out) :: jacobian(size(x, 1), size(parameters)) !! Model Jacobian.
      real(dp) :: decay(size(x, 1)), exponent(size(x, 1)), rate_derivative(size(x, 1))

      exponent = exp(parameters(3))*x(:, 1)**parameters(4)
      decay = exp(-exponent)
      rate_derivative = parameters(2)*decay*exponent
      jacobian(:, 1) = 1.0_dp
      jacobian(:, 2) = -decay
      jacobian(:, 3) = rate_derivative
      jacobian(:, 4) = rate_derivative*log(x(:, 1))
   end subroutine ss_weibull_jacobian

end module r_stats_self_start
