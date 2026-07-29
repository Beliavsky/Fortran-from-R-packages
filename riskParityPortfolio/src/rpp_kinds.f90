! SPDX-License-Identifier: GPL-3.0-only
! Derived from riskParityPortfolio 0.2.2.9000, Copyright Ze Vinicius and Daniel P. Palomar.
module rpp_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter, public :: rpp_huge = huge(1.0_dp) / 16.0_dp
end module rpp_kinds
