! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from OptionPricing 0.1.2 by Wolfgang Hormann and Kemal Dingec.
module optionpricing_types
   use optionpricing_kinds, only : dp
   implicit none
   private

   type, public :: european_result
      real(dp) :: price = 0.0_dp
      real(dp) :: delta = 0.0_dp
      real(dp) :: gamma = 0.0_dp
      real(dp) :: upstream_gamma = 0.0_dp
   end type european_result

   type, public :: greeks_result
      real(dp) :: estimate(3) = 0.0_dp
      real(dp) :: error95(3) = 0.0_dp
      integer :: status = 0
      character(len=160) :: message = ''
   end type greeks_result

   type, public :: moments_result
      real(dp) :: price = 0.0_dp
      real(dp) :: delta = 0.0_dp
      real(dp) :: gamma = 0.0_dp
      integer :: status = 0
      character(len=160) :: message = ''
   end type moments_result

   type, public :: conditional_result
      real(dp) :: y(3) = 0.0_dp
      real(dp) :: controls(6) = 0.0_dp
      integer :: status = 0
      character(len=160) :: message = ''
   end type conditional_result
end module optionpricing_types
