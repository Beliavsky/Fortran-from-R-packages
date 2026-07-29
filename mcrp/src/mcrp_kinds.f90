! SPDX-License-Identifier: GPL-3.0-only
module mcrp_kinds
   implicit none
   private

   integer, parameter, public :: dp = kind(1.0d0)

   integer, parameter, public :: mcrp_success = 0
   integer, parameter, public :: mcrp_invalid_shape = 1
   integer, parameter, public :: mcrp_invalid_argument = 2
   integer, parameter, public :: mcrp_numerical_failure = 3
   integer, parameter, public :: mcrp_max_iterations = 4

end module mcrp_kinds
