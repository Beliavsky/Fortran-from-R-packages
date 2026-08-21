! SPDX-License-Identifier: GPL-2.0-only
! Numerical port of mvtnorm 1.4-2.
module mvtnorm_kinds
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: sqrt_two = sqrt(2.0_dp)
  real(dp), parameter :: log_two_pi = log(2.0_dp*pi)
end module mvtnorm_kinds
