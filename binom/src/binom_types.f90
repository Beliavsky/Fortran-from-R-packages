! SPDX-License-Identifier: GPL-2.0-or-later
module binom_types
   use r_kinds, only: dp
   implicit none
   private
   public :: binom_interval, binom_methods, binom_ok, binom_domain, binom_convergence, binom_bad_method, binom_resource
   integer, parameter :: binom_ok = 0, binom_domain = 1, binom_convergence = 2, binom_bad_method = 3, binom_resource = 4
   character(len=13), parameter :: binom_methods(11) = [character(len=13) :: &
      'agresti-coull', 'asymptotic', 'bayes', 'cloglog', 'exact', 'logit', 'probit', 'profile', 'lrt', 'prop.test', 'wilson']
   type :: binom_interval
      real(dp) :: mean = 0.0_dp
      real(dp) :: lower = 0.0_dp
      real(dp) :: upper = 1.0_dp
      real(dp) :: shape1 = 0.0_dp
      real(dp) :: shape2 = 0.0_dp
      real(dp) :: sig = 0.0_dp
      integer :: status = binom_ok
   end type binom_interval
end module binom_types
