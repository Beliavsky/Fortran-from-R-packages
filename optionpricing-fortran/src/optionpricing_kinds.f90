! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from OptionPricing 0.1.2 by Wolfgang Hormann and Kemal Dingec.
module optionpricing_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter, public :: pi = acos(-1.0_dp)
end module optionpricing_kinds
