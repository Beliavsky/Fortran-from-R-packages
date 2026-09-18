! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the R package wavelets 0.3-0.2 by Eric Aldrich.
! Modern Fortran translation and modifications: 2026-09-10.
module wavelets_kinds
   use, intrinsic :: iso_fortran_env, only : real64
   implicit none
   private

   integer, parameter, public :: dp = real64
end module wavelets_kinds
