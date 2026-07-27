! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_change_point
   use sde_kinds, only : dp
   use sde_interfaces, only : sde_coefficient
   use sde_special, only : normal_pdf
   use sde_nonparametric, only : default_bandwidth
   implicit none
   private

   public :: change_point_result
   public :: detect_change_point

   type :: change_point_result
      integer :: index = 0
      real(dp) :: time = 0.0_dp
      real(dp) :: scale_before = 0.0_dp
      real(dp) :: scale_after = 0.0_dp
      real(dp) :: statistic = 0.0_dp
   end type change_point_result

contains

   subroutine detect_change_point(x, dt, result, t0, drift, diffusion, theta)
      real(dp), intent(in) :: x(:), dt
      type(change_point_result), intent(out) :: result
      real(dp), intent(in), optional :: t0
      procedure(sde_coefficient), optional :: drift, diffusion
      real(dp), intent(in), optional :: theta(:)
      real(dp), allocatable :: z(:), cumulative(:), weights(:)
      real(dp) :: start_time, total, current_stat, max_stat, bw, drift_estimate
      real(dp), allocatable :: local_theta(:)
      integer :: n, len_z, i, best

      n = size(x)
      if (n < 3 .or. dt <= 0.0_dp) error stop "detect_change_point: invalid data or dt"
      if (present(drift) .neqv. present(diffusion)) then
         error stop "detect_change_point: drift and diffusion must be supplied together"
      end if
      start_time = 0.0_dp
      if (present(t0)) start_time = t0
      if (present(theta)) then
         allocate(local_theta(size(theta)))
         local_theta = theta
      else
         allocate(local_theta(0))
      end if
      len_z = n-1
      allocate(z(len_z), cumulative(len_z))
      if (present(drift)) then
         do i = 1, len_z
            z(i) = (x(i+1)-x(i)-drift(0.0_dp, x(i), local_theta)*dt)/ &
               (sqrt(dt)*diffusion(0.0_dp, x(i), local_theta))
         end do
      else
         bw = default_bandwidth(x)
         allocate(weights(len_z))
         do i = 1, len_z
            call kernel_weights_at(x(i), x(1:len_z), bw, weights)
            if (sum(weights) > 0.0_dp) then
               drift_estimate = dot_product(weights, x(2:n)-x(1:len_z))/(dt*sum(weights))
            else
               drift_estimate = 0.0_dp
            end if
            z(i) = (x(i+1)-x(i))/sqrt(dt)-drift_estimate*sqrt(dt)
         end do
      end if

      cumulative(1) = z(1)*z(1)
      do i = 2, len_z
         cumulative(i) = cumulative(i-1)+z(i)*z(i)
      end do
      total = cumulative(len_z)
      if (total <= 0.0_dp) error stop "detect_change_point: residual variance is zero"
      best = 1
      max_stat = -1.0_dp
      do i = 1, len_z-1
         current_stat = abs(real(i, dp)/real(len_z, dp)-cumulative(i)/total)
         if (current_stat > max_stat) then
            max_stat = current_stat
            best = i
         end if
      end do
      result%index = best+1
      result%time = start_time+real(best, dp)*dt
      result%scale_before = sqrt(cumulative(best)/real(best, dp))
      result%scale_after = sqrt((total-cumulative(best))/real(len_z-best, dp))
      result%statistic = max_stat
   end subroutine detect_change_point

   subroutine kernel_weights_at(point, centers, bandwidth, weights)
      real(dp), intent(in) :: point, centers(:), bandwidth
      real(dp), intent(out) :: weights(:)
      integer :: i
      do i = 1, size(centers)
         weights(i) = normal_pdf(point, centers(i), bandwidth)
      end do
   end subroutine kernel_weights_at

end module sde_change_point
