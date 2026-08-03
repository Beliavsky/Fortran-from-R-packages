! SPDX-License-Identifier: GPL-2.0-or-later
module fincal_types
   use fincal_kinds, only : dp
   use fincal_status, only : fincal_ok
   implicit none
   private

   type, public :: inventory_result
      real(dp) :: cost_of_goods = 0.0_dp
      real(dp) :: ending_inventory = 0.0_dp
      integer :: status = fincal_ok
   end type inventory_result

   type, public :: root_result
      real(dp) :: root = 0.0_dp
      real(dp) :: function_value = 0.0_dp
      integer :: iterations = 0
      integer :: status = fincal_ok
   end type root_result
end module fincal_types
