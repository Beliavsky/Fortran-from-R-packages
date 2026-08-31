! SPDX-License-Identifier: GPL-3.0-only
! Derived from computational code in R package pan 2.0.
! Upstream authorship/maintenance: Joseph L. Schafer and Jing Hua Zhao.
! Modern Fortran translation of computational code from R package pan 2.0.
module pan_kinds
   use iso_fortran_env, only : real64, int64
   implicit none
   private

   integer, parameter, public :: dp = real64
   integer, parameter, public :: i8 = int64

end module pan_kinds
