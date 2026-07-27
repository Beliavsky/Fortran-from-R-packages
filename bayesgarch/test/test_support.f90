! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2008-2021 David Ardia
! Modern Fortran translation of computational routines from bayesGARCH.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
module test_support
   use bayesgarch_kinds, only : dp
   implicit none
   private

   public :: assert_close
   public :: assert_true
   public :: test_custom_constraint

   interface assert_close
      module procedure assert_close_scalar
      module procedure assert_close_vector
      module procedure assert_close_matrix
   end interface assert_close

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         write(*, '(a)') "FAILED: " // trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close_scalar(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual
      real(dp), intent(in) :: expected
      real(dp), intent(in) :: tolerance
      character(len=*), intent(in) :: message

      if (abs(actual - expected) > tolerance) then
         write(*, '(a,2(1x,es24.16))') "FAILED: " // trim(message), actual, expected
         error stop 1
      end if
   end subroutine assert_close_scalar

   subroutine assert_close_vector(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual(:)
      real(dp), intent(in) :: expected(:)
      real(dp), intent(in) :: tolerance
      character(len=*), intent(in) :: message

      call assert_true(size(actual) == size(expected), trim(message) // " size")
      if (maxval(abs(actual - expected)) > tolerance) then
         write(*, '(a,1x,es24.16)') "FAILED: " // trim(message) // " max error", &
            maxval(abs(actual - expected))
         error stop 1
      end if
   end subroutine assert_close_vector

   subroutine assert_close_matrix(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual(:, :)
      real(dp), intent(in) :: expected(:, :)
      real(dp), intent(in) :: tolerance
      character(len=*), intent(in) :: message

      call assert_true(all(shape(actual) == shape(expected)), trim(message) // " shape")
      if (maxval(abs(actual - expected)) > tolerance) then
         write(*, '(a,1x,es24.16)') "FAILED: " // trim(message) // " max error", &
            maxval(abs(actual - expected))
         error stop 1
      end if
   end subroutine assert_close_matrix

   logical function test_custom_constraint(psi)
      real(dp), intent(in) :: psi(4)

      test_custom_constraint = psi(1) < 0.06_dp .and. psi(4) < 50.0_dp
   end function test_custom_constraint

end module test_support
