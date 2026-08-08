! SPDX-License-Identifier: GPL-3.0-only
module neighbours_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   integer, parameter, public :: i8 = selected_int_kind(18)
end module neighbours_kinds
