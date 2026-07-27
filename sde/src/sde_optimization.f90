! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_optimization
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use sde_kinds, only : dp
   use sde_interfaces, only : vector_objective
   use sde_utils, only : sort_indices
   implicit none
   private

   public :: optimization_result
   public :: nelder_mead_box

   type :: optimization_result
      real(dp), allocatable :: x(:)
      real(dp) :: value = huge(1.0_dp)
      integer :: iterations = 0
      integer :: evaluations = 0
      integer :: status = 1
      logical :: converged = .false.
   end type optimization_result

contains

   subroutine nelder_mead_box(objective, initial, result, lower, upper, max_iterations, x_tolerance, f_tolerance, initial_step)
      procedure(vector_objective) :: objective
      real(dp), intent(in) :: initial(:)
      type(optimization_result), intent(out) :: result
      real(dp), intent(in), optional :: lower(:)
      real(dp), intent(in), optional :: upper(:)
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: x_tolerance
      real(dp), intent(in), optional :: f_tolerance
      real(dp), intent(in), optional :: initial_step
      real(dp), allocatable :: simplex(:, :), values(:), lo(:), hi(:), centroid(:)
      real(dp), allocatable :: reflected(:), expanded(:), contracted(:), best(:)
      integer, allocatable :: order(:)
      real(dp) :: alpha, gamma, rho, sigma, xtol, ftol, step, trial_value
      real(dp) :: expanded_value, contracted_value
      real(dp) :: simplex_spread, value_spread
      integer :: n, maxiter, iteration, i, best_index, worst_index, second_worst_index

      n = size(initial)
      if (n == 0) error stop "nelder_mead_box: empty parameter vector"
      if (present(lower)) then
         if (size(lower) /= n) error stop "nelder_mead_box: lower bound size mismatch"
      end if
      if (present(upper)) then
         if (size(upper) /= n) error stop "nelder_mead_box: upper bound size mismatch"
      end if

      allocate(simplex(n, n+1), values(n+1), lo(n), hi(n), centroid(n))
      allocate(reflected(n), expanded(n), contracted(n), best(n), order(n+1))
      lo = -huge(1.0_dp)
      hi = huge(1.0_dp)
      if (present(lower)) lo = lower
      if (present(upper)) hi = upper
      if (any(lo > hi)) error stop "nelder_mead_box: lower bound exceeds upper bound"

      maxiter = 2000
      xtol = 1.0e-8_dp
      ftol = 1.0e-10_dp
      step = 0.05_dp
      if (present(max_iterations)) maxiter = max_iterations
      if (present(x_tolerance)) xtol = x_tolerance
      if (present(f_tolerance)) ftol = f_tolerance
      if (present(initial_step)) step = initial_step
      alpha = 1.0_dp
      gamma = 2.0_dp
      rho = 0.5_dp
      sigma = 0.5_dp

      simplex(:, 1) = project(initial)
      do i = 1, n
         simplex(:, i+1) = simplex(:, 1)
         simplex(i, i+1) = simplex(i, i+1)+step*max(1.0_dp, abs(simplex(i, 1)))
         simplex(:, i+1) = project(simplex(:, i+1))
         if (abs(simplex(i, i+1)-simplex(i, 1)) <= &
             epsilon(1.0_dp)*max(1.0_dp, abs(simplex(i, 1)))) then
            simplex(i, i+1) = simplex(i, i+1)-step*max(1.0_dp, abs(simplex(i, 1)))
            simplex(:, i+1) = project(simplex(:, i+1))
         end if
         if (abs(simplex(i, i+1)-simplex(i, 1)) <= &
             epsilon(1.0_dp)*max(1.0_dp, abs(simplex(i, 1))) .and. hi(i) > lo(i)) then
            simplex(i, i+1) = lo(i)+real(i, dp)/real(n+1, dp)*(hi(i)-lo(i))
         end if
      end do

      result%evaluations = 0
      do i = 1, n+1
         values(i) = evaluate(simplex(:, i))
      end do

      do iteration = 1, maxiter
         call sort_indices(values, order)
         best_index = order(1)
         worst_index = order(n+1)
         second_worst_index = order(n)
         best = simplex(:, best_index)
         simplex_spread = 0.0_dp
         do i = 1, n+1
            simplex_spread = max(simplex_spread, maxval(abs(simplex(:, i)-best)))
         end do
         value_spread = maxval(abs(values-values(best_index)))
         if (simplex_spread <= xtol*(1.0_dp+maxval(abs(best))) .and. &
             value_spread <= ftol*(1.0_dp+abs(values(best_index)))) then
            result%converged = .true.
            result%status = 0
            exit
         end if

         centroid = (sum(simplex, dim=2)-simplex(:, worst_index))/real(n, dp)
         reflected = project(centroid+alpha*(centroid-simplex(:, worst_index)))
         trial_value = evaluate(reflected)

         if (trial_value < values(best_index)) then
            expanded = project(centroid+gamma*(reflected-centroid))
            expanded_value = evaluate(expanded)
            if (expanded_value < trial_value) then
               simplex(:, worst_index) = expanded
               values(worst_index) = expanded_value
            else
               simplex(:, worst_index) = reflected
               values(worst_index) = trial_value
            end if
         else if (trial_value < values(second_worst_index)) then
            simplex(:, worst_index) = reflected
            values(worst_index) = trial_value
         else
            if (trial_value < values(worst_index)) then
               contracted = project(centroid+rho*(reflected-centroid))
            else
               contracted = project(centroid-rho*(centroid-simplex(:, worst_index)))
            end if
            contracted_value = evaluate(contracted)
            if (contracted_value < min(trial_value, values(worst_index))) then
               simplex(:, worst_index) = contracted
               values(worst_index) = contracted_value
            else
               do i = 1, n+1
                  if (i == best_index) cycle
                  simplex(:, i) = project(best+sigma*(simplex(:, i)-best))
                  values(i) = evaluate(simplex(:, i))
               end do
            end if
         end if
      end do

      call sort_indices(values, order)
      allocate(result%x(n))
      result%x = simplex(:, order(1))
      result%value = values(order(1))
      result%iterations = min(iteration, maxiter)
      if (.not. result%converged) result%status = 1

   contains

      function project(x) result(y)
         real(dp), intent(in) :: x(:)
         real(dp) :: y(size(x))
         y = max(lo, min(hi, x))
      end function project

      function evaluate(x) result(value)
         real(dp), intent(in) :: x(:)
         real(dp) :: value
         value = objective(x)
         result%evaluations = result%evaluations+1
         if (.not. ieee_is_finite(value)) value = huge(1.0_dp)/16.0_dp
      end function evaluate


   end subroutine nelder_mead_box

end module sde_optimization
