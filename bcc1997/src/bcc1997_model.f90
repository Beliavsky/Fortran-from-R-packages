! SPDX-License-Identifier: GPL-2.0-or-later
! Based on BCC1997 0.1.1, Copyright (C) 2017 Haoran Zhang.
module bcc1997_model
   use bcc1997_kinds, only : dp
   use bcc1997_types, only : bcc_parameters, integration_settings, bcc_result
   use bcc1997_quadrature, only : quadrature_result, integrate_to_infinity
   implicit none
   private

   real(dp), parameter :: pi = acos(-1.0_dp)
   complex(dp), parameter :: imaginary_unit = (0.0_dp, 1.0_dp)

   public :: bcc
   public :: bcc_price
   public :: bcc_price_strikes
   public :: bcc_characteristic_1
   public :: bcc_characteristic_2
   public :: black_scholes_price
   public :: validate_parameters

contains

   function bcc(kappav, kappar, thetav, thetar, sigmav, sigmar, muj, &
      sigmaj, rho, lambda, s0, k, v0, r0, t, settings) result(answer)
      real(dp), intent(in) :: kappav, kappar, thetav, thetar
      real(dp), intent(in) :: sigmav, sigmar, muj, sigmaj, rho, lambda
      real(dp), intent(in) :: s0, k, v0, r0, t
      type(integration_settings), intent(in), optional :: settings
      type(bcc_result) :: answer
      type(bcc_parameters) :: parameters

      parameters = bcc_parameters(kappa_v=kappav, kappa_r=kappar, &
         theta_v=thetav, theta_r=thetar, sigma_v=sigmav, sigma_r=sigmar, &
         mu_j=muj, sigma_j=sigmaj, rho=rho, lambda=lambda, spot=s0, &
         strike=k, variance0=v0, rate0=r0, maturity=t)
      if (present(settings)) then
         answer = bcc_price(parameters, settings)
      else
         answer = bcc_price(parameters)
      end if
   end function bcc

   function bcc_price(parameters, settings) result(answer)
      type(bcc_parameters), intent(in) :: parameters
      type(integration_settings), intent(in), optional :: settings
      type(bcc_result) :: answer
      type(integration_settings) :: local_settings
      type(quadrature_result) :: integral_1, integral_2
      real(dp) :: discount
      logical :: valid
      character(len=160) :: validation_message

      if (present(settings)) then
         local_settings = settings
      else
         local_settings = integration_settings()
      end if

      call validate_parameters(parameters, local_settings, valid, &
         validation_message)
      if (.not. valid) then
         answer%status = 1
         answer%message = validation_message
         return
      end if

      discount = exp(-parameters%rate0 * parameters%maturity)
      if (parameters%maturity <= tiny(1.0_dp)) then
         answer%call = max(parameters%spot - parameters%strike, 0.0_dp)
         answer%put = max(parameters%strike - parameters%spot, 0.0_dp)
         answer%probability1 = merge(1.0_dp, 0.0_dp, &
            parameters%spot > parameters%strike)
         answer%probability2 = answer%probability1
         answer%converged = .true.
         answer%message = 'maturity is zero; intrinsic values returned'
         return
      end if

      integral_1 = integrate_to_infinity(bcc_integrand_1, parameters, local_settings)
      integral_2 = integrate_to_infinity(bcc_integrand_2, parameters, local_settings)

      answer%integral1 = integral_1%value
      answer%integral2 = integral_2%value
      answer%error1 = integral_1%error
      answer%error2 = integral_2%error
      answer%upper_bound1 = integral_1%upper_bound
      answer%upper_bound2 = integral_2%upper_bound
      answer%evaluations1 = integral_1%evaluations
      answer%evaluations2 = integral_2%evaluations
      answer%probability1 = 0.5_dp + integral_1%value / pi
      answer%probability2 = 0.5_dp + integral_2%value / pi
      answer%call = parameters%spot * answer%probability1 - &
         parameters%strike * discount * answer%probability2
      answer%put = -parameters%spot * (1.0_dp - answer%probability1) + &
         parameters%strike * discount * (1.0_dp - answer%probability2)
      answer%converged = integral_1%converged .and. integral_2%converged

      if (answer%converged) then
         answer%message = 'success'
      else
         answer%status = 2
         answer%message = 'improper integration reached its configured limit'
      end if

   end function bcc_price

   function bcc_integrand_1(phi, context) result(value)
      real(dp), intent(in) :: phi
      class(*), intent(in) :: context
      real(dp) :: value
      complex(dp) :: kernel

      select type (parameters => context)
      type is (bcc_parameters)
         kernel = exp(-imaginary_unit * phi * log(parameters%strike)) * &
            bcc_characteristic_1(phi, parameters) / (imaginary_unit * phi)
         value = real(kernel, dp)
      class default
         value = 0.0_dp
      end select
   end function bcc_integrand_1

   function bcc_integrand_2(phi, context) result(value)
      real(dp), intent(in) :: phi
      class(*), intent(in) :: context
      real(dp) :: value
      complex(dp) :: kernel

      select type (parameters => context)
      type is (bcc_parameters)
         kernel = exp(-imaginary_unit * phi * log(parameters%strike)) * &
            bcc_characteristic_2(phi, parameters) / (imaginary_unit * phi)
         value = real(kernel, dp)
      class default
         value = 0.0_dp
      end select
   end function bcc_integrand_2

   subroutine bcc_price_strikes(parameters, strikes, results, settings)
      type(bcc_parameters), intent(in) :: parameters
      real(dp), intent(in) :: strikes(:)
      type(bcc_result), allocatable, intent(out) :: results(:)
      type(integration_settings), intent(in), optional :: settings
      type(bcc_parameters) :: local_parameters
      integer :: i

      allocate(results(size(strikes)))
      local_parameters = parameters
      do i = 1, size(strikes)
         local_parameters%strike = strikes(i)
         if (present(settings)) then
            results(i) = bcc_price(local_parameters, settings)
         else
            results(i) = bcc_price(local_parameters)
         end if
      end do
   end subroutine bcc_price_strikes

   pure function bcc_characteristic_1(phi, parameters) result(value)
      real(dp), intent(in) :: phi
      type(bcc_parameters), intent(in) :: parameters
      complex(dp) :: value
      complex(dp) :: er, ev, qv, one_minus_er, one_minus_ev
      complex(dp) :: a, b, c, d, e, f, g, h
      complex(dp) :: denominator_r, denominator_v, jump_power

      er = sqrt(cmplx(parameters%kappa_r**2, 0.0_dp, dp) - &
         2.0_dp * parameters%sigma_r**2 * imaginary_unit * phi)
      ev = sqrt((parameters%kappa_v - (1.0_dp + imaginary_unit * phi) * &
         parameters%rho * parameters%sigma_v)**2 - imaginary_unit * phi * &
         (1.0_dp + imaginary_unit * phi) * parameters%sigma_v**2)

      one_minus_er = one_minus_exp_minus(er * parameters%maturity)
      one_minus_ev = one_minus_exp_minus(ev * parameters%maturity)
      qv = ev - parameters%kappa_v + (1.0_dp + imaginary_unit * phi) * &
         parameters%rho * parameters%sigma_v

      a = theta_scaled_term(parameters%theta_r, parameters%sigma_r, &
         2.0_dp * log(1.0_dp - (er - parameters%kappa_r) * one_minus_er / &
         (2.0_dp * er)) + (er - parameters%kappa_r) * parameters%maturity)
      b = theta_scaled_term(parameters%theta_v, parameters%sigma_v, &
         2.0_dp * log(1.0_dp - qv * one_minus_ev / (2.0_dp * ev)))
      c = theta_scaled_term(parameters%theta_v, parameters%sigma_v, &
         qv * parameters%maturity)
      d = imaginary_unit * phi * log(parameters%spot)

      denominator_r = 2.0_dp * er - (er - parameters%kappa_r) * one_minus_er
      denominator_v = 2.0_dp * ev - qv * one_minus_ev
      e = 2.0_dp * imaginary_unit * phi * one_minus_er / denominator_r * &
         parameters%rate0

      if (parameters%lambda <= tiny(1.0_dp)) then
         f = (0.0_dp, 0.0_dp)
         g = (0.0_dp, 0.0_dp)
      else
         jump_power = exp(imaginary_unit * phi * log(1.0_dp + &
            parameters%mu_j))
         f = parameters%lambda * (1.0_dp + parameters%mu_j) * &
            parameters%maturity * (jump_power * exp(imaginary_unit * phi * &
            0.5_dp * (1.0_dp + imaginary_unit * phi) * &
            parameters%sigma_j**2) - 1.0_dp)
         g = -parameters%lambda * imaginary_unit * phi * parameters%mu_j * &
            parameters%maturity
      end if

      h = parameters%variance0 * imaginary_unit * phi * &
         (1.0_dp + imaginary_unit * phi) * one_minus_ev / denominator_v
      value = exp(a + b + c + d + e + f + g + h)
   end function bcc_characteristic_1

   pure function bcc_characteristic_2(phi, parameters) result(value)
      real(dp), intent(in) :: phi
      type(bcc_parameters), intent(in) :: parameters
      complex(dp) :: value
      complex(dp) :: er, ev, qv, one_minus_er, one_minus_ev
      complex(dp) :: a, b, c, d, e, f, g, h
      complex(dp) :: denominator_r, denominator_v, jump_power

      er = sqrt(cmplx(parameters%kappa_r**2, 0.0_dp, dp) - &
         2.0_dp * parameters%sigma_r**2 * (imaginary_unit * phi - 1.0_dp))
      ev = sqrt((parameters%kappa_v - imaginary_unit * phi * &
         parameters%rho * parameters%sigma_v)**2 - imaginary_unit * phi * &
         (imaginary_unit * phi - 1.0_dp) * parameters%sigma_v**2)

      one_minus_er = one_minus_exp_minus(er * parameters%maturity)
      one_minus_ev = one_minus_exp_minus(ev * parameters%maturity)
      qv = ev - parameters%kappa_v + imaginary_unit * phi * &
         parameters%rho * parameters%sigma_v

      a = theta_scaled_term(parameters%theta_r, parameters%sigma_r, &
         2.0_dp * log(1.0_dp - (er - parameters%kappa_r) * one_minus_er / &
         (2.0_dp * er)) + (er - parameters%kappa_r) * parameters%maturity)
      b = theta_scaled_term(parameters%theta_v, parameters%sigma_v, &
         2.0_dp * log(1.0_dp - qv * one_minus_ev / (2.0_dp * ev)))
      c = theta_scaled_term(parameters%theta_v, parameters%sigma_v, &
         qv * parameters%maturity)
      d = imaginary_unit * phi * log(parameters%spot) + &
         parameters%rate0 * parameters%maturity

      denominator_r = 2.0_dp * er - (er - parameters%kappa_r) * one_minus_er
      denominator_v = 2.0_dp * ev - qv * one_minus_ev
      e = 2.0_dp * (imaginary_unit * phi - 1.0_dp) * one_minus_er / &
         denominator_r * parameters%rate0

      if (parameters%lambda <= tiny(1.0_dp)) then
         f = (0.0_dp, 0.0_dp)
         g = (0.0_dp, 0.0_dp)
      else
         jump_power = exp(imaginary_unit * phi * log(1.0_dp + &
            parameters%mu_j))
         f = parameters%lambda * parameters%maturity * (jump_power * &
            exp(imaginary_unit * phi * 0.5_dp * &
            (imaginary_unit * phi - 1.0_dp) * parameters%sigma_j**2) - &
            1.0_dp)
         g = -parameters%lambda * imaginary_unit * phi * parameters%mu_j * &
            parameters%maturity
      end if

      h = parameters%variance0 * imaginary_unit * phi * &
         (imaginary_unit * phi - 1.0_dp) * one_minus_ev / denominator_v
      value = exp(a + b + c + d + e + f + g + h)
   end function bcc_characteristic_2

   pure function theta_scaled_term(theta, sigma, expression) result(value)
      real(dp), intent(in) :: theta, sigma
      complex(dp), intent(in) :: expression
      complex(dp) :: value

      if (theta <= tiny(1.0_dp)) then
         value = (0.0_dp, 0.0_dp)
      else
         value = -(theta / sigma**2) * expression
      end if
   end function theta_scaled_term

   pure function one_minus_exp_minus(z) result(value)
      complex(dp), intent(in) :: z
      complex(dp) :: value

      if (abs(z) < 1.0e-5_dp) then
         value = z * (1.0_dp - z * (0.5_dp - z * (1.0_dp / 6.0_dp - &
            z * (1.0_dp / 24.0_dp - z / 120.0_dp))))
      else
         value = 1.0_dp - exp(-z)
      end if
   end function one_minus_exp_minus

   pure subroutine black_scholes_price(spot, strike, rate, volatility, &
      maturity, call, put)
      real(dp), intent(in) :: spot, strike, rate, volatility, maturity
      real(dp), intent(out) :: call, put
      real(dp) :: d1, d2, discount

      if (maturity <= 0.0_dp) then
         call = max(spot - strike, 0.0_dp)
         put = max(strike - spot, 0.0_dp)
         return
      end if

      discount = exp(-rate * maturity)
      if (volatility <= 0.0_dp) then
         call = max(spot - strike * discount, 0.0_dp)
         put = max(strike * discount - spot, 0.0_dp)
         return
      end if

      d1 = (log(spot / strike) + (rate + 0.5_dp * volatility**2) * &
         maturity) / (volatility * sqrt(maturity))
      d2 = d1 - volatility * sqrt(maturity)
      call = spot * normal_cdf(d1) - strike * discount * normal_cdf(d2)
      put = call - spot + strike * discount
   end subroutine black_scholes_price

   pure function normal_cdf(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value

      value = 0.5_dp * erfc(-x / sqrt(2.0_dp))
   end function normal_cdf

   pure subroutine validate_parameters(parameters, settings, valid, message)
      type(bcc_parameters), intent(in) :: parameters
      type(integration_settings), intent(in) :: settings
      logical, intent(out) :: valid
      character(len=*), intent(out) :: message

      valid = .false.
      message = ''

      if (parameters%spot <= 0.0_dp) then
         message = 'spot must be positive'
      else if (parameters%strike <= 0.0_dp) then
         message = 'strike must be positive'
      else if (parameters%maturity < 0.0_dp) then
         message = 'maturity must be nonnegative'
      else if (parameters%variance0 < 0.0_dp) then
         message = 'initial variance must be nonnegative'
      else if (parameters%kappa_v < 0.0_dp .or. parameters%kappa_r < 0.0_dp) then
         message = 'mean-reversion speeds must be nonnegative'
      else if (parameters%theta_v < 0.0_dp .or. parameters%theta_r < 0.0_dp) then
         message = 'long-run variance and rate must be nonnegative'
      else if (parameters%sigma_v <= 0.0_dp .or. parameters%sigma_r <= 0.0_dp) then
         message = 'variance and rate volatilities must be positive'
      else if (parameters%sigma_j < 0.0_dp) then
         message = 'jump volatility must be nonnegative'
      else if (abs(parameters%rho) > 1.0_dp) then
         message = 'rho must be between -1 and 1'
      else if (parameters%lambda < 0.0_dp) then
         message = 'jump intensity must be nonnegative'
      else if (1.0_dp + parameters%mu_j <= 0.0_dp) then
         message = 'one plus the mean jump size must be positive'
      else if (settings%abs_tolerance <= 0.0_dp .or. &
         settings%rel_tolerance <= 0.0_dp) then
         message = 'integration tolerances must be positive'
      else if (settings%panel_width <= 0.0_dp .or. &
         settings%maximum_upper_bound <= settings%panel_width) then
         message = 'invalid integration bounds'
      else if (settings%maximum_depth < 1 .or. settings%tail_panels < 1) then
         message = 'invalid integration iteration limits'
      else
         valid = .true.
         message = 'valid'
      end if
   end subroutine validate_parameters
end module bcc1997_model
