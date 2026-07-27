! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_interfaces
   use sde_kinds, only : dp
   implicit none
   private

   public :: sde_coefficient
   public :: state_function
   public :: vector_objective
   public :: transition_sampler
   public :: moment_function
   public :: estimating_function
   public :: generator_test_function
   public :: martingale_weight_function
   public :: bivariate_function
   public :: scalar_transform

   abstract interface
      pure function sde_coefficient(t, x, theta) result(value)
         import dp
         real(dp), intent(in) :: t
         real(dp), intent(in) :: x
         real(dp), intent(in) :: theta(:)
         real(dp) :: value
      end function sde_coefficient

      pure function state_function(x, theta) result(value)
         import dp
         real(dp), intent(in) :: x
         real(dp), intent(in) :: theta(:)
         real(dp) :: value
      end function state_function

      function vector_objective(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function vector_objective

      function transition_sampler(dt, x0, theta) result(value)
         import dp
         real(dp), intent(in) :: dt
         real(dp), intent(in) :: x0
         real(dp), intent(in) :: theta(:)
         real(dp) :: value
      end function transition_sampler

      subroutine moment_function(y, x, theta, dt, moments)
         import dp
         real(dp), intent(in) :: y
         real(dp), intent(in) :: x
         real(dp), intent(in) :: theta(:)
         real(dp), intent(in) :: dt
         real(dp), intent(out) :: moments(:)
      end subroutine moment_function

      subroutine estimating_function(y, x, theta, values)
         import dp
         real(dp), intent(in) :: y
         real(dp), intent(in) :: x
         real(dp), intent(in) :: theta(:)
         real(dp), intent(out) :: values(:)
      end subroutine estimating_function

      subroutine generator_test_function(index, x, theta, h_x, h_xx)
         import dp
         integer, intent(in) :: index
         real(dp), intent(in) :: x
         real(dp), intent(in) :: theta(:)
         real(dp), intent(out) :: h_x
         real(dp), intent(out) :: h_xx
      end subroutine generator_test_function

      subroutine martingale_weight_function(order, index, x, theta, weight)
         import dp
         integer, intent(in) :: order
         integer, intent(in) :: index
         real(dp), intent(in) :: x
         real(dp), intent(in) :: theta(:)
         real(dp), intent(out) :: weight
      end subroutine martingale_weight_function

      pure function bivariate_function(x, y, theta) result(value)
         import dp
         real(dp), intent(in) :: x
         real(dp), intent(in) :: y
         real(dp), intent(in) :: theta(:)
         real(dp) :: value
      end function bivariate_function

      pure function scalar_transform(x) result(value)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: value
      end function scalar_transform
   end interface

end module sde_interfaces
