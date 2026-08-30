! SPDX-License-Identifier: MIT
! SPDX-FileComment: Shared optional-argument helpers for Fortran translations of R packages.
module r_optional
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private

   integer, parameter :: dp = real64

   public :: optval

   interface optval
      module procedure optval_character, optval_integer, optval_logical, optval_real
   end interface optval

contains

   pure elemental function optval_character(value, default) result(selected)
      !! Returns an optional character value or its fallback.
      character(len=*), intent(in), optional :: value !! Optional value supplied by the caller.
      character(len=*), intent(in) :: default         !! Value returned when `value` is absent.
      character(len=len(default)) :: selected         !! Supplied value or fallback.

      if (present(value)) then
         selected = value
      else
         selected = default
      end if
   end function optval_character

   pure elemental function optval_integer(value, default) result(selected)
      !! Returns an optional integer value or its fallback.
      integer, intent(in), optional :: value !! Optional value supplied by the caller.
      integer, intent(in) :: default         !! Value returned when `value` is absent.
      integer :: selected                    !! Supplied value or fallback.

      if (present(value)) then
         selected = value
      else
         selected = default
      end if
   end function optval_integer

   pure elemental function optval_logical(value, default) result(selected)
      !! Returns an optional logical value or its fallback.
      logical, intent(in), optional :: value !! Optional value supplied by the caller.
      logical, intent(in) :: default         !! Value returned when `value` is absent.
      logical :: selected                    !! Supplied value or fallback.

      if (present(value)) then
         selected = value
      else
         selected = default
      end if
   end function optval_logical

   pure elemental function optval_real(value, default) result(selected)
      !! Returns an optional working-precision real value or its fallback.
      real(dp), intent(in), optional :: value !! Optional value supplied by the caller.
      real(dp), intent(in) :: default         !! Value returned when `value` is absent.
      real(dp) :: selected                    !! Supplied value or fallback.

      if (present(value)) then
         selected = value
      else
         selected = default
      end if
   end function optval_real

end module r_optional
