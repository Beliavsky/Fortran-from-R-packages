! SPDX-License-Identifier: GPL-2.0-or-later
! Real-kind definition for the modern Fortran translation of R package svd 0.5.8.
module svd_kinds
   use iso_fortran_env, only : real64
   implicit none
   private

   integer, parameter, public :: dp = real64
end module svd_kinds
