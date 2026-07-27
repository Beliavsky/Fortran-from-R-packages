! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_information
   use sde_kinds, only : dp, pi
   use sde_interfaces, only : state_function
   use sde_special, only : integrate_adaptive, nan_dp
   use sde_optimization, only : optimization_result, nelder_mead_box
   implicit none
   private

   public :: sde_aic_result
   public :: dc_transition_log_density
   public :: dc_log_likelihood
   public :: sde_aic
   public :: fit_sde_aic

   type :: sde_aic_result
      real(dp), allocatable :: estimate(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      type(optimization_result) :: optimizer
   end type sde_aic_result

contains

   function dc_transition_log_density(x0, x1, dt, theta, drift, diffusion, drift_x, &
         diffusion_x, diffusion_xx) result(value)
      real(dp), intent(in) :: x0, x1, dt, theta(:)
      procedure(state_function) :: drift, diffusion, drift_x, diffusion_x, diffusion_xx
      real(dp) :: value
      real(dp) :: s_transform, h_transform, g_value, sx1

      if (dt <= 0.0_dp) then
         value = nan_dp()
         return
      end if
      sx1 = diffusion(x1, theta)
      if (sx1 <= 0.0_dp) then
         value = nan_dp()
         return
      end if
      s_transform = integrate_adaptive(inverse_sigma, x0, x1, theta)
      h_transform = integrate_adaptive(h_integrand, x0, x1, theta)
      g_value = -0.5_dp*(c1(x0, theta)+c1(x1, theta)+b_transform(x0, theta)*b_transform(x1, theta)/3.0_dp)
      value = -0.5_dp*log(2.0_dp*pi*dt)-log(sx1)-s_transform*s_transform/(2.0_dp*dt)+ &
         h_transform+dt*g_value

   contains

      pure function inverse_sigma(x, local_theta) result(ans)
         real(dp), intent(in) :: x, local_theta(:)
         real(dp) :: ans
         ans = 1.0_dp/diffusion(x, local_theta)
      end function inverse_sigma

      pure function b_transform(x, local_theta) result(ans)
         real(dp), intent(in) :: x, local_theta(:)
         real(dp) :: ans
         ans = drift(x, local_theta)/diffusion(x, local_theta)-0.5_dp*diffusion_x(x, local_theta)
      end function b_transform

      pure function b_transform_x(x, local_theta) result(ans)
         real(dp), intent(in) :: x, local_theta(:)
         real(dp) :: ans, s
         s = diffusion(x, local_theta)
         ans = drift_x(x, local_theta)/s-drift(x, local_theta)*diffusion_x(x, local_theta)/(s*s)- &
            0.5_dp*diffusion_xx(x, local_theta)
      end function b_transform_x

      pure function c1(x, local_theta) result(ans)
         real(dp), intent(in) :: x, local_theta(:)
         real(dp) :: ans, b_value
         b_value = b_transform(x, local_theta)
         ans = b_value*b_value/3.0_dp+0.5_dp*b_transform_x(x, local_theta)*diffusion(x, local_theta)
      end function c1

      pure function h_integrand(x, local_theta) result(ans)
         real(dp), intent(in) :: x, local_theta(:)
         real(dp) :: ans
         ans = b_transform(x, local_theta)/diffusion(x, local_theta)
      end function h_integrand

   end function dc_transition_log_density

   function dc_log_likelihood(x, dt, theta, drift, diffusion, drift_x, diffusion_x, diffusion_xx) result(value)
      real(dp), intent(in) :: x(:), dt, theta(:)
      procedure(state_function) :: drift, diffusion, drift_x, diffusion_x, diffusion_xx
      real(dp) :: value
      integer :: i

      if (size(x) < 2 .or. dt <= 0.0_dp) error stop "dc_log_likelihood: invalid data or dt"
      value = 0.0_dp
      do i = 1, size(x)-1
         value = value+dc_transition_log_density(x(i), x(i+1), dt, theta, drift, diffusion, &
            drift_x, diffusion_x, diffusion_xx)
      end do
   end function dc_log_likelihood

   function sde_aic(x, dt, theta, drift, diffusion, drift_x, diffusion_x, diffusion_xx) result(value)
      real(dp), intent(in) :: x(:), dt, theta(:)
      procedure(state_function) :: drift, diffusion, drift_x, diffusion_x, diffusion_xx
      real(dp) :: value
      value = -2.0_dp*dc_log_likelihood(x, dt, theta, drift, diffusion, drift_x, diffusion_x, diffusion_xx)+ &
         2.0_dp*real(size(theta), dp)
   end function sde_aic

   subroutine fit_sde_aic(x, dt, drift, diffusion, drift_x, diffusion_x, diffusion_xx, initial, result, &
         lower, upper, max_iterations)
      real(dp), intent(in) :: x(:), dt, initial(:)
      procedure(state_function) :: drift, diffusion, drift_x, diffusion_x, diffusion_xx
      type(sde_aic_result), intent(out) :: result
      real(dp), intent(in), optional :: lower(:), upper(:)
      integer, intent(in), optional :: max_iterations
      type(optimization_result) :: opt

      call nelder_mead_box(minimum_contrast, initial, opt, lower=lower, upper=upper, max_iterations=max_iterations)
      result%optimizer = opt
      allocate(result%estimate(size(opt%x)))
      result%estimate = opt%x
      result%log_likelihood = dc_log_likelihood(x, dt, result%estimate, drift, diffusion, drift_x, &
         diffusion_x, diffusion_xx)
      result%aic = -2.0_dp*result%log_likelihood+2.0_dp*real(size(result%estimate), dp)

   contains

      function minimum_contrast(theta) result(value)
         real(dp), intent(in) :: theta(:)
         real(dp) :: value
         real(dp) :: sigma_value, residual, term
         integer :: i
         value = 0.0_dp
         do i = 1, size(x)-1
            sigma_value = diffusion(x(i), theta)
            if (sigma_value <= 0.0_dp) then
               value = huge(1.0_dp)/16.0_dp
               return
            end if
            residual = x(i+1)-x(i)-dt*drift(x(i), theta)
            term = sigma_value+residual*residual/(2.0_dp*dt*sigma_value*sigma_value)
            if (term <= 0.0_dp) then
               value = huge(1.0_dp)/16.0_dp
               return
            end if
            value = value+log(term)
         end do
      end function minimum_contrast

   end subroutine fit_sde_aic

end module sde_information
