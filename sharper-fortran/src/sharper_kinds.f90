! SPDX-License-Identifier: LGPL-3.0-or-later
! Based on SharpeR, copyright 2012-2025 Steven E. Pav.
module sharper_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter, public :: pi = acos(-1.0_dp)
   real(dp), parameter, public :: sqrt_two = sqrt(2.0_dp)
   real(dp), parameter, public :: huge_dp = huge(1.0_dp)
end module sharper_kinds
