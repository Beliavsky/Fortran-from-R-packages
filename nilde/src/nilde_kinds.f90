! SPDX-License-Identifier: GPL-2.0-or-later
module nilde_kinds
   use iso_fortran_env, only : int64, real64
   implicit none
   private
   public :: i8, dp
   integer, parameter :: i8 = int64
   integer, parameter :: dp = real64
end module nilde_kinds
