! SPDX-License-Identifier: GPL-2.0-or-later
module mgcv_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter, public :: pi_dp = acos(-1.0_dp)
end module mgcv_kinds
