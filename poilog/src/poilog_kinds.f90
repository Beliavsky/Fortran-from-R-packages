! SPDX-License-Identifier: GPL-3.0-only
! Derived from the GPL-3 R package poilog by Vidar Grotan and Steinar Engen.
module poilog_kinds
   use, intrinsic :: iso_fortran_env, only : real64
   implicit none
   private
   public :: dp
   integer, parameter :: dp = real64
end module poilog_kinds
