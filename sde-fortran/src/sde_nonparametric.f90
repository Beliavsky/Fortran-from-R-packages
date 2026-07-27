! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_nonparametric
   use sde_kinds, only : dp
   use sde_special, only : normal_pdf
   use sde_utils, only : sample_standard_deviation, linspace
   implicit none
   private

   public :: kernel_estimate
   public :: kernel_drift
   public :: kernel_diffusion
   public :: kernel_density
   public :: default_bandwidth

   type :: kernel_estimate
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: y(:)
      real(dp) :: bandwidth = 0.0_dp
   end type kernel_estimate

contains

   function default_bandwidth(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      if (size(x) < 2) error stop "default_bandwidth: at least two observations are required"
      value = real(size(x), dp)**(-0.2_dp)*sample_standard_deviation(x)
      if (value <= 0.0_dp) value = sqrt(epsilon(1.0_dp))*max(1.0_dp, maxval(abs(x)))
   end function default_bandwidth

   subroutine kernel_drift(x, dt, estimate, bandwidth, n_grid)
      real(dp), intent(in) :: x(:), dt
      type(kernel_estimate), intent(out) :: estimate
      real(dp), intent(in), optional :: bandwidth
      integer, intent(in), optional :: n_grid
      real(dp) :: bw, denominator
      real(dp), allocatable :: weights(:), increments(:)
      integer :: grid_size, i, n

      n = size(x)
      if (n < 2 .or. dt <= 0.0_dp) error stop "kernel_drift: invalid data or dt"
      grid_size = 512
      if (present(n_grid)) grid_size = n_grid
      if (grid_size <= 0) error stop "kernel_drift: n_grid must be positive"
      bw = default_bandwidth(x)
      if (present(bandwidth)) bw = bandwidth
      if (bw <= 0.0_dp) error stop "kernel_drift: bandwidth must be positive"
      estimate%x = linspace(minval(x), maxval(x), grid_size)
      allocate(estimate%y(grid_size), weights(n-1), increments(n-1))
      increments = x(2:n)-x(1:n-1)
      do i = 1, grid_size
         weights = normal_pdf_array(estimate%x(i), x(1:n-1), bw)
         denominator = sum(weights)
         if (denominator > 0.0_dp) then
            estimate%y(i) = dot_product(weights, increments)/(dt*denominator)
         else
            estimate%y(i) = 0.0_dp
         end if
      end do
      estimate%bandwidth = bw
   end subroutine kernel_drift

   subroutine kernel_diffusion(x, dt, estimate, bandwidth, n_grid)
      real(dp), intent(in) :: x(:), dt
      type(kernel_estimate), intent(out) :: estimate
      real(dp), intent(in), optional :: bandwidth
      integer, intent(in), optional :: n_grid
      real(dp) :: bw, denominator, variance_value
      real(dp), allocatable :: weights(:), increments(:)
      integer :: grid_size, i, n

      n = size(x)
      if (n < 2 .or. dt <= 0.0_dp) error stop "kernel_diffusion: invalid data or dt"
      grid_size = 512
      if (present(n_grid)) grid_size = n_grid
      if (grid_size <= 0) error stop "kernel_diffusion: n_grid must be positive"
      bw = default_bandwidth(x)
      if (present(bandwidth)) bw = bandwidth
      if (bw <= 0.0_dp) error stop "kernel_diffusion: bandwidth must be positive"
      estimate%x = linspace(minval(x), maxval(x), grid_size)
      allocate(estimate%y(grid_size), weights(n-1), increments(n-1))
      increments = x(2:n)-x(1:n-1)
      do i = 1, grid_size
         weights = normal_pdf_array(estimate%x(i), x(1:n-1), bw)
         denominator = sum(weights)
         if (denominator > 0.0_dp) then
            variance_value = dot_product(weights, increments*increments)/(dt*denominator)
            estimate%y(i) = sqrt(max(0.0_dp, variance_value))
         else
            estimate%y(i) = 0.0_dp
         end if
      end do
      estimate%bandwidth = bw
   end subroutine kernel_diffusion

   subroutine kernel_density(x, estimate, bandwidth, n_grid, range_min, range_max)
      real(dp), intent(in) :: x(:)
      type(kernel_estimate), intent(out) :: estimate
      real(dp), intent(in), optional :: bandwidth, range_min, range_max
      integer, intent(in), optional :: n_grid
      real(dp) :: bw, lower, upper
      integer :: grid_size, i

      if (size(x) == 0) error stop "kernel_density: empty input"
      grid_size = 512
      if (present(n_grid)) grid_size = n_grid
      bw = default_bandwidth(x)
      if (present(bandwidth)) bw = bandwidth
      if (bw <= 0.0_dp) error stop "kernel_density: bandwidth must be positive"
      lower = minval(x)-3.0_dp*bw
      upper = maxval(x)+3.0_dp*bw
      if (present(range_min)) lower = range_min
      if (present(range_max)) upper = range_max
      estimate%x = linspace(lower, upper, grid_size)
      allocate(estimate%y(grid_size))
      do i = 1, grid_size
         estimate%y(i) = sum(normal_pdf_array(estimate%x(i), x, bw))/real(size(x), dp)
      end do
      estimate%bandwidth = bw
   end subroutine kernel_density

   pure function normal_pdf_array(point, means, sd) result(values)
      real(dp), intent(in) :: point, means(:), sd
      real(dp) :: values(size(means))
      integer :: i
      do i = 1, size(means)
         values(i) = normal_pdf(point, means(i), sd)
      end do
   end function normal_pdf_array

end module sde_nonparametric
