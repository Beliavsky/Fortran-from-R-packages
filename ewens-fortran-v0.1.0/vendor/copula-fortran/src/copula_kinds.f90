! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan
module copula_kinds
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  integer, parameter, public :: i8 = int64
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  real(dp), parameter, public :: sqrt_two = sqrt(2.0_dp)
end module copula_kinds
