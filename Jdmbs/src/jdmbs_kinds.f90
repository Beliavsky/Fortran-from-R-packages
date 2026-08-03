! SPDX-License-Identifier: GPL-2.0-or-later
module jdmbs_kinds
   use, intrinsic :: iso_fortran_env, only : int64
   implicit none
   integer, parameter :: dp = kind(1.0d0)
   public :: dp, int64
end module jdmbs_kinds
