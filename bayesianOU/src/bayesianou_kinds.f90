! SPDX-License-Identifier: MIT
module bayesianou_kinds
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  real(dp), parameter, public :: log_two_pi = log(2.0_dp*pi)
  integer, parameter, public :: status_ok = 0
  integer, parameter, public :: status_bad_input = 1
  integer, parameter, public :: status_not_converged = 2
  integer, parameter, public :: status_singular = 3
  integer, parameter, public :: status_alloc = 4
end module bayesianou_kinds
