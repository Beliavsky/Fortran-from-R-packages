! SPDX-License-Identifier: GPL-2.0-or-later
module moments_status
   implicit none
   private
   integer, parameter, public :: MOMENTS_SUCCESS = 0
   integer, parameter, public :: MOMENTS_INVALID_ARGUMENT = 1
   integer, parameter, public :: MOMENTS_INSUFFICIENT_DATA = 2
   integer, parameter, public :: MOMENTS_DEGENERATE_DATA = 3
   integer, parameter, public :: MOMENTS_NONFINITE_DATA = 4
end module moments_status
