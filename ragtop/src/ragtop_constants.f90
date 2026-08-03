! SPDX-License-Identifier: GPL-2.0-or-later
module ragtop_constants
   implicit none
   private
   integer, parameter, public :: call_option = 1
   integer, parameter, public :: put_option = -1

   integer, parameter, public :: instrument_european_option = 1
   integer, parameter, public :: instrument_american_option = 2
   integer, parameter, public :: instrument_zero_coupon_bond = 3
   integer, parameter, public :: instrument_coupon_bond = 4
   integer, parameter, public :: instrument_callable_bond = 5
   integer, parameter, public :: instrument_convertible_bond = 6

   integer, parameter, public :: ragtop_ok = 0
   integer, parameter, public :: ragtop_invalid_argument = 1
   integer, parameter, public :: ragtop_no_solution = 2
   integer, parameter, public :: ragtop_max_iterations = 3
   integer, parameter, public :: ragtop_singular_system = 4
   integer, parameter, public :: ragtop_allocation_error = 5
end module ragtop_constants
