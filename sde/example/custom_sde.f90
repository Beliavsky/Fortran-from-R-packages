! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
program custom_sde
   use sde
   implicit none

   real(dp), parameter :: theta(3) = [0.8_dp, 0.6_dp, 0.25_dp]
   real(dp), allocatable :: euler_path(:, :), kps_path(:, :)

   call seed_rng(424242_i64)

   call simulate_euler([0.4_dp], 0.0_dp, 0.002_dp, 500, drift, diffusion, theta, &
      euler_path, predictor_corrector=.false.)
   call simulate_kps([0.4_dp], 0.0_dp, 0.002_dp, 500, drift, drift_x, drift_xx, &
      diffusion, diffusion_x, diffusion_xx, theta, kps_path)

   write(*, '(a, f12.6)') "Euler endpoint: ", euler_path(size(euler_path, 1), 1)
   write(*, '(a, f12.6)') "KPS endpoint:   ", kps_path(size(kps_path, 1), 1)

contains

   pure function drift(t, x, local_theta) result(value)
      real(dp), intent(in) :: t, x, local_theta(:)
      real(dp) :: value
      value = local_theta(1)*x-local_theta(2)*x**3+0.0_dp*t
   end function drift

   pure function drift_x(t, x, local_theta) result(value)
      real(dp), intent(in) :: t, x, local_theta(:)
      real(dp) :: value
      value = local_theta(1)-3.0_dp*local_theta(2)*x*x+0.0_dp*t
   end function drift_x

   pure function drift_xx(t, x, local_theta) result(value)
      real(dp), intent(in) :: t, x, local_theta(:)
      real(dp) :: value
      value = -6.0_dp*local_theta(2)*x+0.0_dp*t
   end function drift_xx

   pure function diffusion(t, x, local_theta) result(value)
      real(dp), intent(in) :: t, x, local_theta(:)
      real(dp) :: value
      value = local_theta(3)+0.0_dp*(t+x)
   end function diffusion

   pure function diffusion_x(t, x, local_theta) result(value)
      real(dp), intent(in) :: t, x, local_theta(:)
      real(dp) :: value
      value = 0.0_dp*(t+x+sum(local_theta))
   end function diffusion_x

   pure function diffusion_xx(t, x, local_theta) result(value)
      real(dp), intent(in) :: t, x, local_theta(:)
      real(dp) :: value
      value = 0.0_dp*(t+x+sum(local_theta))
   end function diffusion_xx

end program custom_sde
