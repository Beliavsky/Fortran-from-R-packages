! SPDX-License-Identifier: GPL-3.0-only
! Numerical translation of fHMM 1.4.3.
module fhmm_kinds
   implicit none
   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter, public :: pi = acos(-1.0_dp)
   real(dp), parameter, public :: log_two_pi = log(2.0_dp*pi)
   real(dp), parameter, public :: tiny_prob = 1.0e-300_dp
end module fhmm_kinds
