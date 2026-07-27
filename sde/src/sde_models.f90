! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_models
   use sde_kinds, only : dp
   use sde_special, only : normal_pdf, normal_logpdf, normal_cdf, normal_quantile, &
      lognormal_pdf, lognormal_cdf, lognormal_quantile, gamma_pdf, gamma_logpdf, &
      gamma_cdf, gamma_quantile, safe_expm1, nan_dp
   use sde_distributions, only : noncentral_chi_square_pdf, noncentral_chi_square_logpdf, &
      noncentral_chi_square_cdf, noncentral_chi_square_quantile
   use sde_random, only : random_normal, random_gamma, random_noncentral_chi_square
   implicit none
   private

   public :: ou_conditional_mean
   public :: ou_conditional_variance
   public :: ou_conditional_pdf
   public :: ou_conditional_cdf
   public :: ou_conditional_quantile
   public :: ou_conditional_random
   public :: ou_stationary_pdf
   public :: ou_stationary_cdf
   public :: ou_stationary_quantile
   public :: ou_stationary_random
   public :: gbm_conditional_pdf
   public :: gbm_conditional_cdf
   public :: gbm_conditional_quantile
   public :: gbm_conditional_random
   public :: cir_conditional_pdf
   public :: cir_conditional_cdf
   public :: cir_conditional_quantile
   public :: cir_conditional_random
   public :: cir_stationary_pdf
   public :: cir_stationary_cdf
   public :: cir_stationary_quantile
   public :: cir_stationary_random

contains

   subroutine require_theta(theta, required, name)
      real(dp), intent(in) :: theta(:)
      integer, intent(in) :: required
      character(len=*), intent(in) :: name
      if (size(theta) < required) error stop trim(name)//": theta is too short"
   end subroutine require_theta

   function ou_conditional_mean(dt, x0, theta) result(value)
      real(dp), intent(in) :: dt, x0, theta(:)
      real(dp) :: value
      real(dp) :: level, kappa

      call require_theta(theta, 3, "ou_conditional_mean")
      if (dt < 0.0_dp) error stop "ou_conditional_mean: dt must be nonnegative"
      level = theta(1)
      kappa = theta(2)
      if (abs(kappa*dt) < 1.0e-8_dp) then
         value = x0+level*dt-kappa*x0*dt
      else
         value = level/kappa+(x0-level/kappa)*exp(-kappa*dt)
      end if
   end function ou_conditional_mean

   function ou_conditional_variance(dt, theta) result(value)
      real(dp), intent(in) :: dt, theta(:)
      real(dp) :: value
      real(dp) :: kappa, sigma

      call require_theta(theta, 3, "ou_conditional_variance")
      if (dt < 0.0_dp) error stop "ou_conditional_variance: dt must be nonnegative"
      kappa = theta(2)
      sigma = theta(3)
      if (sigma <= 0.0_dp) error stop "ou_conditional_variance: sigma must be positive"
      if (abs(kappa*dt) < 1.0e-8_dp) then
         value = sigma*sigma*dt
      else
         value = -sigma*sigma*safe_expm1(-2.0_dp*kappa*dt)/(2.0_dp*kappa)
      end if
   end function ou_conditional_variance

   function ou_conditional_pdf(x, dt, x0, theta, log_density) result(value)
      real(dp), intent(in) :: x, dt, x0, theta(:)
      logical, intent(in), optional :: log_density
      real(dp) :: value, mean_value, variance_value
      logical :: log_result

      log_result = .false.
      if (present(log_density)) log_result = log_density
      mean_value = ou_conditional_mean(dt, x0, theta)
      variance_value = ou_conditional_variance(dt, theta)
      if (log_result) then
         value = normal_logpdf(x, mean_value, sqrt(variance_value))
      else
         value = normal_pdf(x, mean_value, sqrt(variance_value))
      end if
   end function ou_conditional_pdf

   function ou_conditional_cdf(x, dt, x0, theta, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: x, dt, x0, theta(:)
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: value
      logical :: lower, log_p

      lower = .true.
      log_p = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_probability)) log_p = log_probability
      value = normal_cdf(x, ou_conditional_mean(dt, x0, theta), &
         sqrt(ou_conditional_variance(dt, theta)), lower)
      if (log_p) value = log(value)
   end function ou_conditional_cdf

   function ou_conditional_quantile(p, dt, x0, theta, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: p, dt, x0, theta(:)
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: value, probability
      logical :: lower, log_p

      lower = .true.
      log_p = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_probability)) log_p = log_probability
      probability = p
      if (log_p) probability = exp(p)
      value = normal_quantile(probability, ou_conditional_mean(dt, x0, theta), &
         sqrt(ou_conditional_variance(dt, theta)), lower)
   end function ou_conditional_quantile

   function ou_conditional_random(dt, x0, theta) result(value)
      real(dp), intent(in) :: dt, x0, theta(:)
      real(dp) :: value
      value = random_normal(ou_conditional_mean(dt, x0, theta), sqrt(ou_conditional_variance(dt, theta)))
   end function ou_conditional_random

   subroutine ou_stationary_parameters(theta, mean_value, sd_value)
      real(dp), intent(in) :: theta(:)
      real(dp), intent(out) :: mean_value, sd_value
      call require_theta(theta, 3, "ou_stationary_parameters")
      if (theta(2) <= 0.0_dp .or. theta(3) <= 0.0_dp) then
         error stop "OU stationary law requires theta(2)>0 and theta(3)>0"
      end if
      mean_value = theta(1)/theta(2)
      sd_value = theta(3)/sqrt(2.0_dp*theta(2))
   end subroutine ou_stationary_parameters

   function ou_stationary_pdf(x, theta, log_density) result(value)
      real(dp), intent(in) :: x, theta(:)
      logical, intent(in), optional :: log_density
      real(dp) :: value, mean_value, sd_value
      logical :: log_result
      call ou_stationary_parameters(theta, mean_value, sd_value)
      log_result = .false.
      if (present(log_density)) log_result = log_density
      if (log_result) then
         value = normal_logpdf(x, mean_value, sd_value)
      else
         value = normal_pdf(x, mean_value, sd_value)
      end if
   end function ou_stationary_pdf

   function ou_stationary_cdf(x, theta, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: x, theta(:)
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: value, mean_value, sd_value
      logical :: lower, log_p
      call ou_stationary_parameters(theta, mean_value, sd_value)
      lower = .true.
      log_p = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_probability)) log_p = log_probability
      value = normal_cdf(x, mean_value, sd_value, lower)
      if (log_p) value = log(value)
   end function ou_stationary_cdf

   function ou_stationary_quantile(p, theta, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: p, theta(:)
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: value, mean_value, sd_value, probability
      logical :: lower, log_p
      call ou_stationary_parameters(theta, mean_value, sd_value)
      lower = .true.
      log_p = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_probability)) log_p = log_probability
      probability = p
      if (log_p) probability = exp(p)
      value = normal_quantile(probability, mean_value, sd_value, lower)
   end function ou_stationary_quantile

   function ou_stationary_random(theta) result(value)
      real(dp), intent(in) :: theta(:)
      real(dp) :: value, mean_value, sd_value
      call ou_stationary_parameters(theta, mean_value, sd_value)
      value = random_normal(mean_value, sd_value)
   end function ou_stationary_random

   subroutine gbm_log_parameters(dt, x0, theta, meanlog, sdlog)
      real(dp), intent(in) :: dt, x0, theta(:)
      real(dp), intent(out) :: meanlog, sdlog
      call require_theta(theta, 2, "gbm_log_parameters")
      if (dt < 0.0_dp .or. x0 <= 0.0_dp .or. theta(2) <= 0.0_dp) then
         error stop "GBM law requires dt>=0, x0>0, and sigma>0"
      end if
      meanlog = log(x0)+(theta(1)-0.5_dp*theta(2)**2)*dt
      sdlog = theta(2)*sqrt(dt)
   end subroutine gbm_log_parameters

   function gbm_conditional_pdf(x, dt, x0, theta, log_density) result(value)
      real(dp), intent(in) :: x, dt, x0, theta(:)
      logical, intent(in), optional :: log_density
      real(dp) :: value, meanlog, sdlog
      logical :: log_result
      call gbm_log_parameters(dt, x0, theta, meanlog, sdlog)
      log_result = .false.
      if (present(log_density)) log_result = log_density
      if (log_result) then
         if (x <= 0.0_dp) then
            value = -huge(1.0_dp)
         else
            value = normal_logpdf(log(x), meanlog, sdlog)-log(x)
         end if
      else
         value = lognormal_pdf(x, meanlog, sdlog)
      end if
   end function gbm_conditional_pdf

   function gbm_conditional_cdf(x, dt, x0, theta, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: x, dt, x0, theta(:)
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: value, meanlog, sdlog
      logical :: lower, log_p
      call gbm_log_parameters(dt, x0, theta, meanlog, sdlog)
      lower = .true.
      log_p = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_probability)) log_p = log_probability
      value = lognormal_cdf(x, meanlog, sdlog, lower)
      if (log_p) value = log(value)
   end function gbm_conditional_cdf

   function gbm_conditional_quantile(p, dt, x0, theta, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: p, dt, x0, theta(:)
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: value, meanlog, sdlog, probability
      logical :: lower, log_p
      call gbm_log_parameters(dt, x0, theta, meanlog, sdlog)
      lower = .true.
      log_p = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_probability)) log_p = log_probability
      probability = p
      if (log_p) probability = exp(p)
      value = lognormal_quantile(probability, meanlog, sdlog, lower)
   end function gbm_conditional_quantile

   function gbm_conditional_random(dt, x0, theta) result(value)
      real(dp), intent(in) :: dt, x0, theta(:)
      real(dp) :: value, meanlog, sdlog
      call gbm_log_parameters(dt, x0, theta, meanlog, sdlog)
      value = exp(random_normal(meanlog, sdlog))
   end function gbm_conditional_random

   subroutine cir_conditional_parameters(dt, x0, theta, scale_c, df, ncp)
      real(dp), intent(in) :: dt, x0, theta(:)
      real(dp), intent(out) :: scale_c, df, ncp
      real(dp) :: decay
      call require_theta(theta, 3, "cir_conditional_parameters")
      if (dt <= 0.0_dp .or. x0 < 0.0_dp .or. any(theta(1:3) <= 0.0_dp)) then
         error stop "CIR law requires dt>0, x0>=0, and positive parameters"
      end if
      decay = exp(-theta(2)*dt)
      scale_c = 2.0_dp*theta(2)/((1.0_dp-decay)*theta(3)**2)
      ncp = 2.0_dp*scale_c*x0*decay
      df = 4.0_dp*theta(1)/theta(3)**2
   end subroutine cir_conditional_parameters

   function cir_conditional_pdf(x, dt, x0, theta, log_density) result(value)
      real(dp), intent(in) :: x, dt, x0, theta(:)
      logical, intent(in), optional :: log_density
      real(dp) :: value, scale_c, df, ncp, transformed
      logical :: log_result
      call cir_conditional_parameters(dt, x0, theta, scale_c, df, ncp)
      log_result = .false.
      if (present(log_density)) log_result = log_density
      transformed = 2.0_dp*scale_c*x
      if (log_result) then
         value = log(2.0_dp*scale_c)+noncentral_chi_square_logpdf(transformed, df, ncp)
      else
         value = 2.0_dp*scale_c*noncentral_chi_square_pdf(transformed, df, ncp)
      end if
   end function cir_conditional_pdf

   function cir_conditional_cdf(x, dt, x0, theta, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: x, dt, x0, theta(:)
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: value, scale_c, df, ncp
      logical :: lower, log_p
      call cir_conditional_parameters(dt, x0, theta, scale_c, df, ncp)
      lower = .true.
      log_p = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_probability)) log_p = log_probability
      value = noncentral_chi_square_cdf(2.0_dp*scale_c*x, df, ncp, lower)
      if (log_p) value = log(value)
   end function cir_conditional_cdf

   function cir_conditional_quantile(p, dt, x0, theta, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: p, dt, x0, theta(:)
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: value, scale_c, df, ncp, probability
      logical :: lower, log_p
      call cir_conditional_parameters(dt, x0, theta, scale_c, df, ncp)
      lower = .true.
      log_p = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_probability)) log_p = log_probability
      probability = p
      if (log_p) probability = exp(p)
      value = noncentral_chi_square_quantile(probability, df, ncp, lower)/(2.0_dp*scale_c)
   end function cir_conditional_quantile

   function cir_conditional_random(dt, x0, theta) result(value)
      real(dp), intent(in) :: dt, x0, theta(:)
      real(dp) :: value, scale_c, df, ncp
      call cir_conditional_parameters(dt, x0, theta, scale_c, df, ncp)
      value = random_noncentral_chi_square(df, ncp)/(2.0_dp*scale_c)
   end function cir_conditional_random

   subroutine cir_stationary_parameters(theta, shape, scale)
      real(dp), intent(in) :: theta(:)
      real(dp), intent(out) :: shape, scale
      call require_theta(theta, 3, "cir_stationary_parameters")
      if (any(theta(1:3) <= 0.0_dp)) error stop "CIR stationary law requires positive parameters"
      shape = 2.0_dp*theta(1)/theta(3)**2
      scale = theta(3)**2/(2.0_dp*theta(2))
   end subroutine cir_stationary_parameters

   function cir_stationary_pdf(x, theta, log_density) result(value)
      real(dp), intent(in) :: x, theta(:)
      logical, intent(in), optional :: log_density
      real(dp) :: value, shape, scale
      logical :: log_result
      call cir_stationary_parameters(theta, shape, scale)
      log_result = .false.
      if (present(log_density)) log_result = log_density
      if (log_result) then
         value = gamma_logpdf(x, shape, scale)
      else
         value = gamma_pdf(x, shape, scale)
      end if
   end function cir_stationary_pdf

   function cir_stationary_cdf(x, theta, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: x, theta(:)
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: value, shape, scale
      logical :: lower, log_p
      call cir_stationary_parameters(theta, shape, scale)
      lower = .true.
      log_p = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_probability)) log_p = log_probability
      value = gamma_cdf(x, shape, scale, lower)
      if (log_p) value = log(value)
   end function cir_stationary_cdf

   function cir_stationary_quantile(p, theta, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: p, theta(:)
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: value, shape, scale, probability
      logical :: lower, log_p
      call cir_stationary_parameters(theta, shape, scale)
      lower = .true.
      log_p = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_probability)) log_p = log_probability
      probability = p
      if (log_p) probability = exp(p)
      value = gamma_quantile(probability, shape, scale, lower)
   end function cir_stationary_quantile

   function cir_stationary_random(theta) result(value)
      real(dp), intent(in) :: theta(:)
      real(dp) :: value, shape, scale
      call cir_stationary_parameters(theta, shape, scale)
      value = random_gamma(shape, scale)
   end function cir_stationary_random

end module sde_models
