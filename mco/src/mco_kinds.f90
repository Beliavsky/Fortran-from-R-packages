! SPDX-License-Identifier: GPL-2.0-only
module mco_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter, public :: pi = acos(-1.0_dp)
end module mco_kinds
