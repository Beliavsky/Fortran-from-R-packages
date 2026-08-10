! SPDX-License-Identifier: GPL-3.0-only
! Translation of computational code from the R package scs (GPL-3),
! which bundles SCS by Brendan O'Donoghue under the MIT license.
module scs_kinds
   use, intrinsic :: iso_fortran_env, only : real64, int32
   implicit none
   private
   public :: dp, i4
   integer, parameter :: dp = real64
   integer, parameter :: i4 = int32
end module scs_kinds
