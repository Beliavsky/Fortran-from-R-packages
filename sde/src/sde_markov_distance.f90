! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_markov_distance
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use sde_kinds, only : dp
   use sde_utils, only : interpolate_missing
   implicit none
   private

   public :: markov_operator_distance
   public :: bspline_basis_matrix

contains

   subroutine markov_operator_distance(x, n_basis, distances, operators, range_min, range_max, spline_order)
      real(dp), intent(in) :: x(:, :)
      integer, intent(in) :: n_basis
      real(dp), allocatable, intent(out) :: distances(:, :)
      real(dp), allocatable, intent(out), optional :: operators(:, :, :)
      real(dp), intent(in), optional :: range_min, range_max
      integer, intent(in), optional :: spline_order
      real(dp), allocatable :: clean(:, :), basis(:, :), local_operators(:, :, :)
      real(dp) :: lower, upper
      integer :: n, n_series, order, series, i, j, t

      n = size(x, 1)
      n_series = size(x, 2)
      if (n < 2 .or. n_series < 1) error stop "markov_operator_distance: invalid matrix"
      order = 10
      if (present(spline_order)) order = spline_order
      if (order < 1 .or. n_basis < order) then
         error stop "markov_operator_distance: n_basis must be at least spline_order"
      end if
      lower = minval(x, mask=.not. ieee_is_nan(x))
      upper = maxval(x, mask=.not. ieee_is_nan(x))
      if (present(range_min)) lower = range_min
      if (present(range_max)) upper = range_max
      if (upper <= lower) error stop "markov_operator_distance: degenerate range"

      allocate(clean(n, n_series), local_operators(n_basis, n_basis, n_series))
      clean = x
      do series = 1, n_series
         if (any(ieee_is_nan(clean(:, series)))) call interpolate_missing(clean(:, series))
         call bspline_basis_matrix(clean(:, series), lower, upper, n_basis, order, basis)
         local_operators(:, :, series) = 0.0_dp
         do i = 1, n_basis
            do j = 1, n_basis
               do t = 1, n-1
                  local_operators(i, j, series) = local_operators(i, j, series)+ &
                     basis(t, i)*basis(t+1, j)+basis(t, j)*basis(t+1, i)
               end do
               local_operators(i, j, series) = local_operators(i, j, series)/(2.0_dp*real(n, dp))
            end do
         end do
      end do

      allocate(distances(n_series, n_series))
      distances = 0.0_dp
      do i = 1, n_series-1
         do j = i+1, n_series
            distances(i, j) = sum(abs(local_operators(:, :, i)-local_operators(:, :, j)))
            distances(j, i) = distances(i, j)
         end do
      end do
      if (present(operators)) then
         allocate(operators(n_basis, n_basis, n_series))
         operators = local_operators
      end if
   end subroutine markov_operator_distance

   subroutine bspline_basis_matrix(x, lower, upper, n_basis, order, basis)
      real(dp), intent(in) :: x(:), lower, upper
      integer, intent(in) :: n_basis, order
      real(dp), allocatable, intent(out) :: basis(:, :)
      real(dp), allocatable :: knots(:), current(:), next(:)
      real(dp) :: denominator_left, denominator_right
      integer :: n_knots, n_interior, obs, i, level

      if (upper <= lower .or. order < 1 .or. n_basis < order) then
         error stop "bspline_basis_matrix: invalid basis specification"
      end if
      n_knots = n_basis+order
      n_interior = n_basis-order
      allocate(knots(n_knots), current(n_basis+order), next(n_basis+order))
      knots(1:order) = lower
      do i = 1, n_interior
         knots(order+i) = lower+(upper-lower)*real(i, dp)/real(n_interior+1, dp)
      end do
      knots(n_basis+1:n_knots) = upper
      allocate(basis(size(x), n_basis))
      basis = 0.0_dp

      do obs = 1, size(x)
         current = 0.0_dp
         do i = 1, n_basis+order-1
            if ((x(obs) >= knots(i) .and. x(obs) < knots(i+1)) .or. &
                (abs(x(obs)-upper) <= 4.0_dp*spacing(upper) .and. i == n_basis)) then
               current(i) = 1.0_dp
            end if
         end do
         do level = 2, order
            next = 0.0_dp
            do i = 1, n_basis+order-level
               denominator_left = knots(i+level-1)-knots(i)
               denominator_right = knots(i+level)-knots(i+1)
               if (denominator_left > 0.0_dp) then
                  next(i) = next(i)+(x(obs)-knots(i))/denominator_left*current(i)
               end if
               if (denominator_right > 0.0_dp) then
                  next(i) = next(i)+(knots(i+level)-x(obs))/denominator_right*current(i+1)
               end if
            end do
            current = next
         end do
         basis(obs, :) = current(1:n_basis)
      end do
   end subroutine bspline_basis_matrix

end module sde_markov_distance
