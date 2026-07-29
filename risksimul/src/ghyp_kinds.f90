! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Reused under GPL-2/GPL-3 from the ghyp-fortran numerical implementation.
! Derived from ghyp 1.6.5 by Marc Weibel, David Luethi, and Henriette-Elise Breymann.
module ghyp_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   integer, parameter, public :: i8 = selected_int_kind(18)
   real(dp), parameter, public :: pi = acos(-1.0_dp)
   real(dp), parameter, public :: sqrt_two_pi = sqrt(2.0_dp*pi)
   real(dp), parameter, public :: log_two_pi = log(2.0_dp*pi)
end module ghyp_kinds
