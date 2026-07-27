! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_density
   use sde_kinds, only : dp, pi
   use sde_interfaces, only : sde_coefficient
   use sde_special, only : normal_pdf, normal_logpdf, safe_log1p, expm1_over_x, safe_expm1, nan_dp
   implicit none
   private

   public :: transition_density_euler
   public :: transition_density_elerian
   public :: transition_density_kessler
   public :: transition_density_ozaki
   public :: transition_density_shoji
   public :: ozaki_moments
   public :: shoji_moments

contains

   function transition_density_euler(x, t, x0, t0, theta, drift, diffusion, log_density) result(value)
      real(dp), intent(in) :: x, t, x0, t0, theta(:)
      procedure(sde_coefficient) :: drift, diffusion
      logical, intent(in), optional :: log_density
      real(dp) :: value, dt, mean_value, sd_value
      logical :: log_result

      dt = t-t0
      if (dt <= 0.0_dp) then
         value = nan_dp()
         return
      end if
      sd_value = sqrt(dt)*abs(diffusion(t0, x0, theta))
      mean_value = x0+drift(t0, x0, theta)*dt
      log_result = .false.
      if (present(log_density)) log_result = log_density
      if (log_result) then
         value = normal_logpdf(x, mean_value, sd_value)
      else
         value = normal_pdf(x, mean_value, sd_value)
      end if
   end function transition_density_euler

   function transition_density_elerian(x, t, x0, t0, theta, drift, diffusion, diffusion_x, &
         log_density) result(value)
      real(dp), intent(in) :: x, t, x0, t0, theta(:)
      procedure(sde_coefficient) :: drift, diffusion, diffusion_x
      logical, intent(in), optional :: log_density
      real(dp) :: value, dt, a, b, z, c, root, log_cosh, s, sx, log_value
      logical :: log_result

      dt = t-t0
      if (dt <= 0.0_dp) then
         value = nan_dp()
         return
      end if
      s = diffusion(t0, x0, theta)
      sx = diffusion_x(t0, x0, theta)
      a = 0.5_dp*s*sx*dt
      if (abs(a) <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(x0)) .or. abs(sx) <= tiny(1.0_dp)) then
         value = transition_density_euler(x, t, x0, t0, theta, drift, diffusion, log_density)
         return
      end if
      b = -s/(2.0_dp*sx)+x0+drift(t0, x0, theta)*dt-a
      z = (x-b)/a
      if (z <= 0.0_dp .or. abs(s) <= tiny(1.0_dp)) then
         if (present(log_density)) then
            if (log_density) then
               value = -huge(1.0_dp)
               return
            end if
         end if
         value = 0.0_dp
         return
      end if
      c = 1.0_dp/(s*s*dt)
      root = sqrt(c*z)
      log_cosh = abs(root)+safe_log1p(exp(-2.0_dp*abs(root)))-log(2.0_dp)
      log_value = -0.5_dp*(c+z)+log_cosh-0.5_dp*log(z)-log(abs(a))-log(2.0_dp*pi)
      log_result = .false.
      if (present(log_density)) log_result = log_density
      if (log_result) then
         value = log_value
      else
         value = exp(log_value)
      end if
   end function transition_density_elerian

   function transition_density_kessler(x, t, x0, t0, theta, drift, drift_x, drift_xx, &
         diffusion, diffusion_x, diffusion_xx, log_density) result(value)
      real(dp), intent(in) :: x, t, x0, t0, theta(:)
      procedure(sde_coefficient) :: drift, drift_x, drift_xx
      procedure(sde_coefficient) :: diffusion, diffusion_x, diffusion_xx
      logical, intent(in), optional :: log_density
      real(dp) :: value, dt, mu, mu1, mu2, sg, sg1, sg2, mean_value, variance_value
      logical :: log_result

      dt = t-t0
      if (dt <= 0.0_dp) then
         value = nan_dp()
         return
      end if
      mu = drift(t0, x0, theta)
      mu1 = drift_x(t0, x0, theta)
      mu2 = drift_xx(t0, x0, theta)
      sg = diffusion(t0, x0, theta)
      sg1 = diffusion_x(t0, x0, theta)
      sg2 = diffusion_xx(t0, x0, theta)
      mean_value = x0+mu*dt+(mu*mu1+0.5_dp*sg*sg*mu2)*dt*dt/2.0_dp
      variance_value = x0*x0+(2.0_dp*mu*x0+sg*sg)*dt+ &
         (2.0_dp*mu*(mu1*x0+mu+sg*sg1)+sg*sg*(mu2*x0+2.0_dp*mu1+sg1*sg1+sg*sg2))* &
         dt*dt/2.0_dp-mean_value*mean_value
      if (variance_value <= 0.0_dp) then
         value = nan_dp()
         return
      end if
      log_result = .false.
      if (present(log_density)) log_result = log_density
      if (log_result) then
         value = normal_logpdf(x, mean_value, sqrt(variance_value))
      else
         value = normal_pdf(x, mean_value, sqrt(variance_value))
      end if
   end function transition_density_kessler

   subroutine ozaki_moments(t0, x0, dt, theta, drift, drift_x, diffusion, mean_value, variance_value)
      real(dp), intent(in) :: t0, x0, dt, theta(:)
      procedure(sde_coefficient) :: drift, drift_x, diffusion
      real(dp), intent(out) :: mean_value, variance_value
      real(dp) :: d, l, s, k, adjustment

      d = drift(t0, x0, theta)
      l = drift_x(t0, x0, theta)
      s = diffusion(t0, x0, theta)
      mean_value = x0+d*dt*expm1_over_x(l*dt)
      if (abs(x0) > sqrt(epsilon(1.0_dp))) then
         adjustment = d*dt*expm1_over_x(l*dt)/x0
         if (adjustment > -1.0_dp) then
            k = safe_log1p(adjustment)/dt
         else
            k = l
         end if
      else
         k = l
      end if
      variance_value = s*s*dt*expm1_over_x(2.0_dp*k*dt)
   end subroutine ozaki_moments

   function transition_density_ozaki(x, t, x0, t0, theta, drift, drift_x, diffusion, &
         log_density) result(value)
      real(dp), intent(in) :: x, t, x0, t0, theta(:)
      procedure(sde_coefficient) :: drift, drift_x, diffusion
      logical, intent(in), optional :: log_density
      real(dp) :: value, mean_value, variance_value
      logical :: log_result

      if (t <= t0) then
         value = nan_dp()
         return
      end if
      call ozaki_moments(t0, x0, t-t0, theta, drift, drift_x, diffusion, mean_value, variance_value)
      if (variance_value <= 0.0_dp) then
         value = nan_dp()
         return
      end if
      log_result = .false.
      if (present(log_density)) log_result = log_density
      if (log_result) then
         value = normal_logpdf(x, mean_value, sqrt(variance_value))
      else
         value = normal_pdf(x, mean_value, sqrt(variance_value))
      end if
   end function transition_density_ozaki

   subroutine shoji_moments(t0, x0, dt, theta, drift, drift_x, drift_xx, drift_t, diffusion, &
         mean_value, variance_value)
      real(dp), intent(in) :: t0, x0, dt, theta(:)
      procedure(sde_coefficient) :: drift, drift_x, drift_xx, drift_t, diffusion
      real(dp), intent(out) :: mean_value, variance_value
      real(dp) :: d, l, s, m, z

      d = drift(t0, x0, theta)
      l = drift_x(t0, x0, theta)
      s = diffusion(t0, x0, theta)
      m = 0.5_dp*s*s*drift_xx(t0, x0, theta)+drift_t(t0, x0, theta)
      z = l*dt
      mean_value = x0+d*dt*expm1_over_x(z)+m*dt*dt*expm1_minus_x_over_x2(z)
      variance_value = s*s*dt*expm1_over_x(2.0_dp*z)
   end subroutine shoji_moments

   function transition_density_shoji(x, t, x0, t0, theta, drift, drift_x, drift_xx, drift_t, &
         diffusion, log_density) result(value)
      real(dp), intent(in) :: x, t, x0, t0, theta(:)
      procedure(sde_coefficient) :: drift, drift_x, drift_xx, drift_t, diffusion
      logical, intent(in), optional :: log_density
      real(dp) :: value, mean_value, variance_value
      logical :: log_result

      if (t <= t0) then
         value = nan_dp()
         return
      end if
      call shoji_moments(t0, x0, t-t0, theta, drift, drift_x, drift_xx, drift_t, diffusion, &
         mean_value, variance_value)
      if (variance_value <= 0.0_dp) then
         value = nan_dp()
         return
      end if
      log_result = .false.
      if (present(log_density)) log_result = log_density
      if (log_result) then
         value = normal_logpdf(x, mean_value, sqrt(variance_value))
      else
         value = normal_pdf(x, mean_value, sqrt(variance_value))
      end if
   end function transition_density_shoji

   pure function expm1_minus_x_over_x2(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      if (abs(x) < 1.0e-5_dp) then
         value = 0.5_dp+x/6.0_dp+x*x/24.0_dp+x*x*x/120.0_dp
      else
         value = (safe_expm1(x)-x)/(x*x)
      end if
   end function expm1_minus_x_over_x2

end module sde_density
