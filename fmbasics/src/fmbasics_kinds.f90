! SPDX-License-Identifier: GPL-2.0-only
module fmbasics_kinds
   implicit none
   private

   integer, parameter, public :: dp = kind(1.0d0)
   integer, parameter, public :: FM_OK = 0
   integer, parameter, public :: FM_INVALID_ARGUMENT = 1
   integer, parameter, public :: FM_SIZE_MISMATCH = 2
   integer, parameter, public :: FM_DOMAIN_ERROR = 3
   integer, parameter, public :: FM_NOT_CONVERGED = 4
   integer, parameter, public :: FM_IO_ERROR = 5
   integer, parameter, public :: FM_UNSUPPORTED = 6
end module fmbasics_kinds
