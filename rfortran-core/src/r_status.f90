! SPDX-License-Identifier: MIT
module r_status
   implicit none
   private

   integer, parameter, public :: r_ok = 0
   integer, parameter, public :: r_invalid_input = 1
   integer, parameter, public :: r_no_data = 2
   integer, parameter, public :: r_singular = 3

end module r_status
