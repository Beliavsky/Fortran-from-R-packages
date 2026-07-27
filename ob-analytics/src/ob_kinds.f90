! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of obAnalytics.
! Copyright (C) 2015,2016 Philip Stubbings.
module ob_kinds
   use, intrinsic :: iso_fortran_env, only : int64, real64
   implicit none
   private
   integer, parameter, public :: dp = real64
   integer, parameter, public :: i8 = int64
end module ob_kinds
