! SPDX-License-Identifier: GPL-2.0-or-later
module highs_kinds
   use, intrinsic :: iso_c_binding, only : c_int, c_int64_t, c_double
   implicit none
   private
   integer, parameter, public :: dp = c_double
   integer, parameter, public :: highs_int = c_int
   integer, parameter, public :: highs_int64 = c_int64_t
end module highs_kinds
