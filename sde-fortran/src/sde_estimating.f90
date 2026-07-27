! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_estimating
   use sde_kinds, only : dp
   use sde_interfaces, only : estimating_function, sde_coefficient, generator_test_function, &
      state_function, martingale_weight_function
   use sde_optimization, only : optimization_result, nelder_mead_box
   implicit none
   private

   public :: estimating_result
   public :: fit_simple_estimating
   public :: fit_generator_estimating
   public :: fit_linear_martingale
   public :: evaluate_simple_estimating
   public :: evaluate_generator_estimating
   public :: evaluate_linear_martingale

   type :: estimating_result
      real(dp), allocatable :: estimate(:)
      real(dp), allocatable :: equations(:)
      real(dp) :: objective = huge(1.0_dp)
      type(optimization_result) :: optimizer
   end type estimating_result

contains

   subroutine evaluate_simple_estimating(x, theta, n_equations, estimating, equations)
      real(dp), intent(in) :: x(:), theta(:)
      integer, intent(in) :: n_equations
      procedure(estimating_function) :: estimating
      real(dp), intent(out) :: equations(:)
      real(dp), allocatable :: values(:)
      integer :: i

      if (size(x) < 2 .or. size(equations) /= n_equations) then
         error stop "evaluate_simple_estimating: invalid dimensions"
      end if
      allocate(values(n_equations))
      equations = 0.0_dp
      do i = 1, size(x)-1
         call estimating(x(i+1), x(i), theta, values)
         equations = equations+values
      end do
   end subroutine evaluate_simple_estimating

   subroutine fit_simple_estimating(x, n_equations, estimating, initial, result, lower, upper, max_iterations)
      real(dp), intent(in) :: x(:), initial(:)
      integer, intent(in) :: n_equations
      procedure(estimating_function) :: estimating
      type(estimating_result), intent(out) :: result
      real(dp), intent(in), optional :: lower(:), upper(:)
      integer, intent(in), optional :: max_iterations
      type(optimization_result) :: opt

      call nelder_mead_box(objective, initial, opt, lower=lower, upper=upper, max_iterations=max_iterations)
      result%optimizer = opt
      allocate(result%estimate(size(opt%x)), result%equations(n_equations))
      result%estimate = opt%x
      call evaluate_simple_estimating(x, result%estimate, n_equations, estimating, result%equations)
      result%objective = sum(result%equations**2)

   contains

      function objective(theta) result(value)
         real(dp), intent(in) :: theta(:)
         real(dp) :: value
         real(dp) :: equations_local(n_equations)
         call evaluate_simple_estimating(x, theta, n_equations, estimating, equations_local)
         value = sum(equations_local**2)
      end function objective

   end subroutine fit_simple_estimating

   subroutine evaluate_generator_estimating(x, theta, n_functions, drift, diffusion, h_derivatives, equations)
      real(dp), intent(in) :: x(:), theta(:)
      integer, intent(in) :: n_functions
      procedure(sde_coefficient) :: drift, diffusion
      procedure(generator_test_function) :: h_derivatives
      real(dp), intent(out) :: equations(:)
      real(dp) :: hx, hxx, d, s
      integer :: i, j

      if (size(x) == 0 .or. size(equations) /= n_functions) then
         error stop "evaluate_generator_estimating: invalid dimensions"
      end if
      equations = 0.0_dp
      do i = 1, size(x)
         d = drift(0.0_dp, x(i), theta)
         s = diffusion(0.0_dp, x(i), theta)
         do j = 1, n_functions
            call h_derivatives(j, x(i), theta, hx, hxx)
            equations(j) = equations(j)+d*hx+0.5_dp*s*s*hxx
         end do
      end do
   end subroutine evaluate_generator_estimating

   subroutine fit_generator_estimating(x, n_functions, drift, diffusion, h_derivatives, initial, result, &
         lower, upper, max_iterations)
      real(dp), intent(in) :: x(:), initial(:)
      integer, intent(in) :: n_functions
      procedure(sde_coefficient) :: drift, diffusion
      procedure(generator_test_function) :: h_derivatives
      type(estimating_result), intent(out) :: result
      real(dp), intent(in), optional :: lower(:), upper(:)
      integer, intent(in), optional :: max_iterations
      type(optimization_result) :: opt

      call nelder_mead_box(objective, initial, opt, lower=lower, upper=upper, max_iterations=max_iterations)
      result%optimizer = opt
      allocate(result%estimate(size(opt%x)), result%equations(n_functions))
      result%estimate = opt%x
      call evaluate_generator_estimating(x, result%estimate, n_functions, drift, diffusion, &
         h_derivatives, result%equations)
      result%objective = sum(result%equations**2)

   contains

      function objective(theta) result(value)
         real(dp), intent(in) :: theta(:)
         real(dp) :: value
         real(dp) :: equations_local(n_functions)
         call evaluate_generator_estimating(x, theta, n_functions, drift, diffusion, &
            h_derivatives, equations_local)
         value = sum(equations_local**2)
      end function objective

   end subroutine fit_generator_estimating

   subroutine evaluate_linear_martingale(x, theta, conditional_mean, conditional_variance, weights, &
         equations, use_second_order)
      real(dp), intent(in) :: x(:), theta(:)
      procedure(state_function) :: conditional_mean, conditional_variance
      procedure(martingale_weight_function) :: weights
      real(dp), intent(out) :: equations(:)
      logical, intent(in), optional :: use_second_order
      real(dp) :: h1, h2, weight
      logical :: second_order
      integer :: i, j, n_parameters

      n_parameters = size(equations)
      if (size(x) < 2 .or. n_parameters == 0) then
         error stop "evaluate_linear_martingale: invalid dimensions"
      end if
      second_order = .false.
      if (present(use_second_order)) second_order = use_second_order
      equations = 0.0_dp
      do i = 1, size(x)-1
         h1 = x(i+1)-conditional_mean(x(i), theta)
         h2 = h1*h1-conditional_variance(x(i), theta)
         do j = 1, n_parameters
            call weights(1, j, x(i), theta, weight)
            equations(j) = equations(j)+weight*h1
            if (second_order) then
               call weights(2, j, x(i), theta, weight)
               equations(j) = equations(j)+weight*h2
            end if
         end do
      end do
   end subroutine evaluate_linear_martingale

   subroutine fit_linear_martingale(x, conditional_mean, conditional_variance, weights, initial, result, &
         lower, upper, use_second_order, max_iterations)
      real(dp), intent(in) :: x(:), initial(:)
      procedure(state_function) :: conditional_mean, conditional_variance
      procedure(martingale_weight_function) :: weights
      type(estimating_result), intent(out) :: result
      real(dp), intent(in), optional :: lower(:), upper(:)
      logical, intent(in), optional :: use_second_order
      integer, intent(in), optional :: max_iterations
      type(optimization_result) :: opt
      logical :: second_order
      integer :: n_parameters

      n_parameters = size(initial)
      second_order = .false.
      if (present(use_second_order)) second_order = use_second_order
      call nelder_mead_box(objective, initial, opt, lower=lower, upper=upper, max_iterations=max_iterations)
      result%optimizer = opt
      allocate(result%estimate(n_parameters), result%equations(n_parameters))
      result%estimate = opt%x
      call evaluate_linear_martingale(x, result%estimate, conditional_mean, conditional_variance, &
         weights, result%equations, second_order)
      result%objective = sum(result%equations**2)

   contains

      function objective(theta) result(value)
         real(dp), intent(in) :: theta(:)
         real(dp) :: value
         real(dp) :: equations_local(n_parameters)
         call evaluate_linear_martingale(x, theta, conditional_mean, conditional_variance, &
            weights, equations_local, second_order)
         value = sum(equations_local**2)
      end function objective

   end subroutine fit_linear_martingale

end module sde_estimating
